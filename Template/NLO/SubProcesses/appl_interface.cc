#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <pineappl.h>

#include "orders.h"

/*
  fNLO mode of aMCatNLO
*/

// Declare grids
std::vector<pineappl_grid*> grid_obs;

// translates an index from the range [0, __amp_split_size) to the indices need by `fill_grid`
std::vector<std::vector<int>> translation_tables;

// Maximum number of (pairwise) suprocesses
const int __max_nproc__ = 121;

// Information defined at the generation (configuration) step, that does
// not vary event by event
extern "C" struct
{
    int amp_split_size; // Maximum number of coupling-combinations
    int qcdpower[__amp_split_size]; // Power of alpha_s for each amp_split
    int qedpower[__amp_split_size]; // Power of alpha for each amp_split
} appl_common_fixed_;

// Map of the PDF combinations from aMC@NLO - structure for each
// "subprocess" i, has some number nproc[i] pairs of parton
// combinations. To be used together with the info in appl_flavmap.
extern "C" struct
{
    int lumimap[__max_nproc__][__max_nproc__][2]; // (paired) subprocesses per combination
    int nproc[__max_nproc__]; // number of separate (pairwise) subprocesses for this combination
    int nlumi; // overall number of combinations ( 0 < nlumi <= __mxpdflumi__ )
} appl_common_lumi_;

// Event weights, kinematics, etc. that are different event by event
extern "C" struct
{
    double x1[4], x2[4];
    double muF2[4], muR2[4], muQES2[4];
    double W0[__amp_split_size][4], WR[__amp_split_size][4];
    double WF[__amp_split_size][4], WB[__amp_split_size][4];
    int flavmap[4];
} appl_common_weights_;

// Parameters of the grids.
// These parameters can optionally be singularly specified by the user,
// but if no specification is given, the code will use the default values.
extern "C" struct
{
    double Q2min, Q2max;
    double xmin, xmax;
    int nQ2, Q2order;
    int nx, xorder;
} appl_common_grid_;

// Parameters of the histograms
extern "C" struct
{
    double www_histo, norm_histo;
    double obs_histo;
    double obs_min, obs_max;
    double obs_bins[101];
    int obs_nbins;
    int itype_histo;
    int amp_pos;
    int obs_num;
} appl_common_histokin_;

// Event weight and cross section
extern "C" struct
{
    double event_weight, vegaswgt;
    double xsec12, xsec11, xsec20;
} appl_common_reco_;

// Check if a file exists
bool file_exists(const std::string& s)
{
    if (std::FILE* testfile = std::fopen(s.c_str(), "r"))
    {
        std::fclose(testfile);
        return true;
    }
    else
    {
        return false;
    }
}

// Banner
std::string Banner()
{
    return "\n"
        "    █████╗ ███╗   ███╗ ██████╗██████╗ ██╗      █████╗ ███████╗████████╗\n"
        "   ██╔══██╗████╗ ████║██╔════╝██╔══██╗██║     ██╔══██╗██╔════╝╚══██╔══╝\n"
        "   ███████║██╔████╔██║██║     ██████╔╝██║     ███████║███████╗   ██║\n"
        "   ██╔══██║██║╚██╔╝██║██║     ██╔══██╗██║     ██╔══██║╚════██║   ██║\n"
        "   ██║  ██║██║ ╚═╝ ██║╚██████╗██████╔╝███████╗██║  ██║███████║   ██║\n"
        "   ╚═╝  ╚═╝╚═╝     ╚═╝ ╚═════╝╚═════╝ ╚══════╝╚═╝  ╚═╝╚══════╝   ╚═╝\n";
}

extern "C" void appl_init_()
{
    // Grid Initialization and definition of the observables.
    // Construct the input file name according to its position in the
    // vector "grid_obs".
    std::size_t const index = grid_obs.size();
    std::string const grid_filename_in = "grid_obs_" + std::to_string(index) + "_in.root";

    // Check that the grid file exists. If so read the grid from the file,
    // otherwise create a new grid from scratch.
    if (file_exists(grid_filename_in))
    {
        // Open the existing grid
        grid_obs.emplace_back(pineappl_grid_read(grid_filename_in.c_str()));

        std::vector<int> subgrid_params(4 * pineappl_grid_get_subgrids(grid_obs.back()));
        pineappl_grid_get_subgrid_params(grid_obs.back(), subgrid_params.data());

        translation_tables.emplace_back();
        auto& translation_table = translation_tables.back();

        // when loading the grids, the stored orders might be sorted differently
        for (int i = 0; i != appl_common_fixed_.amp_split_size; ++i)
        {
            std::size_t index = 0;
            bool found = false;

            while (!found)
            {
                auto const alphas = subgrid_params.at(4 * index + 0);
                auto const alpha  = subgrid_params.at(4 * index + 1);
                auto const logxir = subgrid_params.at(4 * index + 2);
                auto const logxif = subgrid_params.at(4 * index + 3);

                if ((alphas == appl_common_fixed_.qcdpower[i] / 2) &&
                    (alpha  == appl_common_fixed_.qedpower[i] / 2) &&
                    (logxir == 0) &&
                    (logxif == 0))
                {
                    found = true;
                }
                else
                {
                    ++index;
                }
            }

            // TODO: can this fail?
            assert( found );

            translation_table.push_back(index);
        }
    }
    // If the grid does not exist, book it after having defined all the
    // relevant parameters.
    else
    {
        int lo_power = 9999;
        int nlo_power = 0;

        for (int i = 0; i != appl_common_fixed_.amp_split_size; ++i)
        {
            int sum = appl_common_fixed_.qcdpower[i] + appl_common_fixed_.qedpower[i];

            lo_power = std::min(lo_power, sum);
            nlo_power = std::max(nlo_power, sum);
        }

        // TODO: are there any situations is which there are NLOs but no LOs?

        // we assume that there is always at least one LO, and zero or one NLO
        assert((nlo_power == (lo_power + 2)) || (nlo_power == lo_power));

        std::vector<int> subgrid_params;

        translation_tables.emplace_back();
        translation_tables.back().reserve(appl_common_fixed_.amp_split_size);

        for (int i = 0; i != appl_common_fixed_.amp_split_size; ++i)
        {
            int const qcd = appl_common_fixed_.qcdpower[i];
            int const qed = appl_common_fixed_.qedpower[i];
            int const sum = qcd + qed;

            translation_tables.back().push_back(subgrid_params.size() / 4);

            if (sum == lo_power)
            {
                // WB
                subgrid_params.insert(subgrid_params.end(), { qcd / 2, qed / 2, 0, 0 });
            }
            else if (sum == nlo_power)
            {
                // W0
                subgrid_params.insert(subgrid_params.end(), { qcd / 2, qed / 2, 0, 0 });
                // WR
                subgrid_params.insert(subgrid_params.end(), { qcd / 2, qed / 2, 1, 0 });
                // WF
                subgrid_params.insert(subgrid_params.end(), { qcd / 2, qed / 2, 0, 1 });
            }
        }

        // Define the settings for the interpolation in x and Q2.
        // These are common to all the grids computed.
        // If values larger than zero (i.e. set by the user) are found the default
        // settings are replaced with the new ones.
        auto* storage = pineappl_storage_new("papplgrid_f2");

        if (appl_common_grid_.nQ2 > 0)
        {
            pineappl_storage_set_int(storage, "nq2", appl_common_grid_.nQ2);
        }
        // Max and min value of Q2
        if (appl_common_grid_.Q2min > 0.0)
        {
            pineappl_storage_set_double(storage, "q2min", appl_common_grid_.Q2min);
        }
        if (appl_common_grid_.Q2max > 0.0)
        {
            pineappl_storage_set_double(storage, "q2max", appl_common_grid_.Q2max);
        }
        // Order of the polynomial interpolation in Q2
        if (appl_common_grid_.Q2order > 0)
        {
            pineappl_storage_set_int(storage, "q2order", appl_common_grid_.Q2order);
        }
        // Number of points for the x interpolation
        if (appl_common_grid_.nx > 0)
        {
            pineappl_storage_set_int(storage, "nx", appl_common_grid_.nx);
        }
        // Min and max value of x
        if (appl_common_grid_.xmin > 0.0)
        {
            pineappl_storage_set_double(storage, "xmin", appl_common_grid_.xmin);
        }
        if (appl_common_grid_.xmax > 0.0)
        {
            pineappl_storage_set_double(storage, "xmax", appl_common_grid_.xmax);
        }
        // Order of the polynomial interpolation in x
        if (appl_common_grid_.xorder > 0)
        {
            pineappl_storage_set_int(storage, "xorder", appl_common_grid_.xorder);
        }

        // Set up the PDF luminosities
        auto* lumi = pineappl_lumi_new();

        // Loop over parton luminosities
        for (int ilumi = 0; ilumi < appl_common_lumi_.nlumi; ilumi++)
        {
            int nproc = appl_common_lumi_.nproc[ilumi];

            std::vector<int> pdg_ids;
            pdg_ids.reserve(2 * nproc);

            for (int iproc = 0; iproc != nproc; ++iproc)
            {
                pdg_ids.push_back(appl_common_lumi_.lumimap[ilumi][iproc][0]);
                pdg_ids.push_back(appl_common_lumi_.lumimap[ilumi][iproc][1]);
            }

            pineappl_lumi_add(lumi, nproc, pdg_ids.data(), nullptr);
        }

        // Use the reweighting function
        pineappl_storage_set_bool(storage, "reweight", true);
        // Add documentation
        pineappl_storage_set_string(storage, "documentation", Banner().c_str());

        // Create a grid with the binning given in the "obsbins[Nbins+1]" array
        grid_obs.push_back(pineappl_grid_new(
            lumi,
            pineappl_subgrid_format::as_a_logxir_logxif,
            subgrid_params.size() / 4,
            subgrid_params.data(),
            appl_common_histokin_.obs_nbins,
            appl_common_histokin_.obs_bins,
            storage
        ));

        pineappl_storage_delete(storage);
        pineappl_lumi_delete(lumi);
    }
}

extern "C" void appl_fill_()
{
    // Check that itype ranges from 1 to 5.
    int itype = appl_common_histokin_.itype_histo;
    if ((itype < 1) || (itype > 5))
    {
        std::cout << "amcblast ERROR: Invalid value of itype = " << itype << std::endl;
        std::exit(-10);
    }

    // this is the second index of the WB/R/F/0 arrays
    int index = appl_common_histokin_.amp_pos - 1;

    // aMC@NLO weights. Four grids, ordered as {W0,WR,WF,WB}.
    double(&W0)[4] = appl_common_weights_.W0[index];
    double(&WR)[4] = appl_common_weights_.WR[index];
    double(&WF)[4] = appl_common_weights_.WF[index];
    double(&WB)[4] = appl_common_weights_.WB[index];

    int ilumi;
    int nlumi = appl_common_lumi_.nlumi;
    double ttol = 1e-100;
    double x1, x2;
    double scale2;
    double obs = appl_common_histokin_.obs_histo;
    // Weight vector whose size is the total number of subprocesses
    std::vector<double> weight(nlumi, 0.0);

    // Histogram number
    int nh = appl_common_histokin_.obs_num - 1;

    // translate (index,nh) -> index of the APPLgrid
    int const grid_index = translation_tables.at(nh).at(index);

    int k;

    switch (itype)
    {
    case 1:
        k = 0;
        break;

    case 2:
    case 3:
        k = 1;
        break;

    case 4:
        k = 2;
        break;

    case 5:
        k = 3;
        break;

    default:
        assert( false );
    }

    x1 = appl_common_weights_.x1[k];
    x2 = appl_common_weights_.x2[k];

    static std::vector<std::vector<std::vector<double>>> x1Saved(5, std::vector<std::vector<double>>(
        grid_obs.size(), std::vector<double>(pineappl_grid_get_subgrids(grid_obs[nh]), 0.0)));
    static std::vector<std::vector<std::vector<double>>> x2Saved(5, std::vector<std::vector<double>>(
        grid_obs.size(), std::vector<double>(pineappl_grid_get_subgrids(grid_obs[nh]), 0.0)));

    if (x1 == x1Saved[itype - 1][nh][grid_index] && x2 == x2Saved[itype - 1][nh][grid_index])
    {
        return;
    }
    else
    {
        x1Saved[itype - 1][nh][grid_index] = x1;
        x2Saved[itype - 1][nh][grid_index] = x2;
    }

    scale2 = appl_common_weights_.muF2[k];
    ilumi = appl_common_weights_.flavmap[k] - 1;

    if (x1 < 0.0 || x1 > 1.0 || x2 < 0.0 || x2 > 1.0)
    {
        std::cout << "amcblast ERROR: Invalid value of x1 and/or x2 = " << x1 << " " << x2
                  << std::endl;
        std::exit(-10);
    }

    if (x1 == 0.0 && x2 == 0.0)
    {
        return;
    }

    if (std::fabs(W0[k]) > ttol)
    {
        // W0
        weight.at(ilumi) = W0[k];
        pineappl_grid_fill(grid_obs[nh], x1, x2, scale2, obs, &weight[0], grid_index + 0);
        weight.at(ilumi) = 0.0;
    }

    if (std::fabs(WR[k]) > ttol)
    {
        // WR
        weight.at(ilumi) = WR[k];
        pineappl_grid_fill(grid_obs[nh], x1, x2, scale2, obs, &weight[0], grid_index + 1);
        weight.at(ilumi) = 0.0;
    }

    if (std::fabs(WF[k]) > ttol)
    {
        // WF
        weight.at(ilumi) = WF[k];
        pineappl_grid_fill(grid_obs[nh], x1, x2, scale2, obs, &weight[0], grid_index + 2);
        weight.at(ilumi) = 0.0;
    }

    if (std::fabs(WB[k]) > ttol)
    {
        // WB
        weight.at(ilumi) = WB[k];
        pineappl_grid_fill(grid_obs[nh], x1, x2, scale2, obs, &weight[0], grid_index);
        weight.at(ilumi) = 0.0;
    }
}

extern "C" void appl_term_()
{
    // Histogram number
    int const nh = appl_common_histokin_.obs_num - 1;

    // Normalize the grid by conversion factor and number of runs
    pineappl_grid_scale(grid_obs[nh], 389379660.0 * appl_common_histokin_.norm_histo);

    // Write grid to file
    pineappl_grid_write(grid_obs[nh], ("grid_obs_" + std::to_string(nh) + "_out.root").c_str());

    pineappl_grid_delete(grid_obs[nh]);
}
