!
!     Module      : DiscreteSampler
!     Author      : Valentin Hirschi
!     Date        : 29.10.2014
!     Destriction : 
!              A relatively simple and flexible module to do
!              sampling of discrete dimensions for Monte-Carlo
!              purposes.
!
!     List of public subroutines and usage :
!
!
!     DS_register_dimension(name, n_bins)
!       ::  Register a new dimension with its name and number of bins
!
!     DS_remove_dimension(name)
!       ::  Removes and clear the dimension of input name
!
!     DS_print_global_info(name|index|void)
!       ::  Print global information on a registered information, using
!       ::  either its name or index. or all if none is selected
!     
!     DS_clear
!       ::  Reinitializes the module and all grid data
!
!     DS_binID(integer_ID)
!       ::  Creates an object of binID type from an integer ID. Notice
!       ::  that you can also use the assignment operator directly.
!
!     DS_add_bin(dimension_name, (binID|integerID|void))
!       ::  Add one bin to dimension of name dimension_name.
!       ::  The user can add an ID to the added bin or just use the
!       ::  default sequential labeling
!    
!     DS_remove_bin(dimension_name, (binID|integerID))
!       ::  Remove a bin in a dimension by either specfiying its index
!       ::  in this dimension or its binID (remember you can create a
!       ::  binID on the flight from an integer using the function
!       ::  DS_binID )
!
!     DS_add_entry(dimension_name, (binID|integerID), weight)
!       ::  Add a new weight to a certan bin (characterized by either
!       ::  its binID or the integer of the binID)
!
!     DS_update_grid((dim_name|void))
!       ::  Update the reference grid of the dimension dim_name or 
!       ::  update all of them at the same time without argument.
!       ::  It uses the running grid for this and it reinitilizes it.
!
!     DS_write_grid((file_name|stream_id), (dim_name|void))
!       :: Append to file 'file_name' or a given stream the data for 
!       :: the current reference grid for dimension dim_name or
!       :: or all of them if called without dim_name.
!       :: Notice that
!
!     DS_load_grid((file_name|stream_id), (dim_name|void))
!       :: Reset the run_grid for dimension dim_name and loads in
!       :: the data obtained from file 'file_name' or stream_id for this
!       :: dimension or all of them if called without dim_name.
!       :: The data is loaded in the running grid only (which is
!       :: re-initialized before this), so you have to call 
!       :: 'DS_update_grid' if you want it pased to the ref_grid.
!
!       --- DONE UP TO HERE ---
!
!     DS_get_point(dim_name, random_variable, 
!                       (binIDPicked|integerIDPicked), jacobian_weight)
!       :: From a given random variable in [0.0,1.0] and a dimension
!       :: name, this subroutine returns the picked bin or index and
!       :: the Jacobian weight associated.
!       :: Jacobian = nbins_in_dim * normalized_wgt_in_selected_bin
!       ::          = wgt_in_selected_bin / wgt_averaged_over_all_bins
!
      module DiscreteSampler

      use StringCast

      logical    DS_verbose
      parameter (DS_verbose=.FALSE.)

!     This parameter sets how large must be the sampling bar when
!     displaying information about a dimension
      integer samplingBarWidth
      parameter (samplingBarWidth=80)

!     Attributes identifying a bin
!     For now just an integer
      type binID
        integer id
      endtype
!     And an easy way to create a binIDs
      interface assignment (=)
        module procedure  binID_from_binID
        module procedure  binID_from_integer
      end interface assignment (=)
!     Define and easy way of comparing binIDs 
      interface operator (==)
        module procedure  equal_binID
      end interface operator (==)

!     Information relevant to a bin
      type bin
        real*8 weight
!       Sum of the squared weights, to compute the variance
        real*8 weight_sqr
!       Sum of the absolute value of the weights put in this bin,
!       necessary when in presence of negative weights.
        real*8 abs_weight
        integer n_entries
!       Practical to be able to identify a bin by its id
        type(binID) bid
      endtype

!     Define and easy way of adding Bins 
      interface operator (+)
        module procedure  DS_combine_two_bins
      end interface operator (+)

      type sampledDimension
!       These are the reference weights, for the grid currently used
!       and controlling the sampling
        type(bin) , dimension(:), allocatable    :: bins
!       Keep track of the norm (i.e. sum of all weights) and the total
!       number of points for ease and optimisation purpose
        real*8                                   :: norm
!       The sum of the absolute value of the weight in each bin
        real*8                                   :: abs_norm
!       The sum of the squared weights in each bin
        real*8                                   :: norm_sqr     
        integer                                  :: n_tot_entries
!       A handy way of referring to the dimension by its name rather than
!       an index.
        character, dimension(:), allocatable     :: dimension_name
      endtype sampledDimension

!     This stores the overall discrete reference grid
      type(sampledDimension), dimension(:), allocatable :: ref_grid

!       This is the running grid, whose weights are being updated for each point
!       but not yet used for the sampling. The user must call the 'update'
!       function for the running grid to be merged to the reference one.
      type(sampledDimension), dimension(:), allocatable :: run_grid

      interface DS_add_entry
        module procedure DS_add_entry_with_BinID
        module procedure DS_add_entry_with_BinIntID
      end interface DS_add_entry ! DS_add_entry

      interface DS_print_global_info
        module procedure DS_print_dim_global_info_from_name
        module procedure DS_print_dim_global_info_from_void
      end interface ! DS_print_dim_global_info

      interface DS_add_bin
        module procedure DS_add_bin_with_binID
        module procedure DS_add_bin_with_IntegerID
        module procedure DS_add_bin_with_void
      end interface ! DS_add_bin

      interface DS_remove_bin
        module procedure DS_remove_bin_withIntegerID
        module procedure DS_remove_bin_withBinID
      end interface ! DS_remove_bin

      interface DS_get_bin
        module procedure DS_get_bin_from_binID
        module procedure DS_get_bin_from_binID_and_dimName
      end interface ! DS_get_bin

      interface DS_update_grid
        module procedure DS_update_grid_with_dim_name
        module procedure DS_update_all_grids
        module procedure DS_update_grid_with_dim_index
      end interface ! DS_update_grid

      interface DS_write_grid
        module procedure DS_write_grid_with_filename
        module procedure DS_write_grid_with_streamID
      end interface ! DS_write_grid

      interface DS_load_grid
        module procedure DS_load_grid_with_filename
        module procedure DS_load_grid_with_streamID
      end interface ! DS_load_grid

      interface DS_dim_index
        module procedure DS_dim_index_default
        module procedure DS_dim_index_with_force
        module procedure DS_dim_index_default_with_chararray
        module procedure DS_dim_index_with_force_with_chararray
      end interface ! DS_dim_index

      interface DS_bin_index
        module procedure DS_bin_index_default
        module procedure DS_bin_index_with_force
      end interface ! DS_bin_index


      contains

!       ---------------------------------------------------------------
!       This subroutine is simply the logger of this module
!       ---------------------------------------------------------------

        subroutine DS_Logger(msg)
        implicit none
!         
!         Subroutine arguments
!         
          character(len=*), intent(in)        :: msg

          if (DS_verbose) write(*,*) msg

        end subroutine DS_Logger

!       ---------------------------------------------------------------
!       This subroutine clears the module and reinitialize all data 
!       ---------------------------------------------------------------
        subroutine DS_clear()
          call DS_deallocate_grid(ref_grid)
          call DS_deallocate_grid(run_grid)
        end subroutine DS_clear

        subroutine DS_deallocate_grid(grid)
          integer i
          type(sampledDimension), dimension(:), allocatable,
     &                                            intent(inout) :: grid
          if (allocated(grid)) then
            do i = 1,size(grid)
              if (allocated(grid(i)%bins)) then
                deallocate(grid(i)%bins)
              endif
              if (allocated(grid(i)%dimension_name)) then
                deallocate(grid(i)%dimension_name)
              endif
            enddo
            deallocate(grid)            
          endif
        end subroutine DS_deallocate_grid

!       ---------------------------------------------------------------
!       This subroutine takes care of registering a new dimension in
!       the DSampler module by characterizin it by its name and number
!       of bins.
!       ---------------------------------------------------------------
        subroutine DS_register_dimension(dim_name,n_bins,all_grids)
        implicit none
!         
!         Subroutine arguments
!        
          integer , intent(in)                :: n_bins
          character(len=*), intent(in)        :: dim_name
          logical , optional                  :: all_grids
!
!         Local variables
!
          logical                             :: do_all_grids
!
!         Begin code
!
          if (present(all_grids)) then
            do_all_grids = all_grids
          else
            do_all_grids = .True.
          endif
          if (do_all_grids) then
            call DS_add_dimension_to_grid(ref_grid, dim_name, n_bins)
          endif
          call DS_add_dimension_to_grid(run_grid, dim_name, n_bins)

          call DS_Logger("DiscreteSampler:: Successfully registered "//
     $    "dimension '"//dim_name//"' ("//TRIM(toStr(n_bins))//' bins)')

        end subroutine DS_register_dimension

!       ---------------------------------------------------------------
!       This subroutine registers a dimension to a given grid 
!       ---------------------------------------------------------------
        subroutine DS_add_dimension_to_grid(grid, dim_name, n_bins)
        implicit none
!         
!         Subroutine arguments
!
          type(sampledDimension), dimension(:), allocatable,
     &      intent(inout)                          :: grid
          integer , intent(in)                     :: n_bins
          character(len=*), intent(in)             :: dim_name
!
!         Local variables
!
          integer                                           :: dim_index
          type(sampledDimension), dimension(:), allocatable :: tmp
          integer i
!
!         Begin code
!
          if(allocated(grid)) then
            dim_index = DS_dim_index(grid,dim_name,.True.)
            if (dim_index.ne.-1) then
               write(*,*) 'DiscreteSampler:: Error, the dimension'//
     $              " with name '"//dim_name//"' is already registered."
               stop 1
            endif
          endif

!         Either allocate the discrete grids or append a dimension 
          if (allocated(grid)) then
            allocate(tmp(size(grid)))
            do i=1, size(grid)
              call DS_copy_dimension(grid(i), tmp(i))
            enddo
            call DS_deallocate_grid(grid)
            allocate(grid(size(tmp)+1))
            do i=1, size(tmp)
              call DS_copy_dimension(tmp(i), grid(i))
            enddo
            call DS_deallocate_grid(tmp)
          else
            allocate(grid(1))
          endif
!         Now we can fill in the appended element with the
!         characteristics of the dimension which must be added
          allocate(grid(size(grid))%bins(n_bins))
          allocate(grid(size(grid))%dimension_name(len(dim_name)))
!         Initialize the values of the grid with default
          call DS_initialize_dimension(grid(size(grid)))
!         For the string assignation, I have to it character by
!         character.
          do i=1, len(dim_name)
            grid(size(grid))%dimension_name(i) = dim_name(i:i)
          enddo

        end subroutine DS_add_dimension_to_grid

!       ----------------------------------------------------------------------
!       Copy a dimension from source to target, making sure to allocate  
!       ----------------------------------------------------------------------
        subroutine DS_copy_dimension(source, trget)
          type(sampledDimension), intent(out)   :: trget
          type(sampledDimension), intent(in)    :: source
          integer i

          if (allocated(trget%bins)) then
            deallocate(trget%bins)
          endif
          allocate(trget%bins(size(source%bins)))
          do i=1,size(source%bins)
            call DS_copy_bin(source%bins(i),trget%bins(i))
          enddo
          if (allocated(trget%dimension_name)) then
            deallocate(trget%dimension_name)
          endif
          allocate(trget%dimension_name(size(source%dimension_name)))
          do i=1, size(source%dimension_name)
            trget%dimension_name(i) = source%dimension_name(i)
          enddo
          trget%norm          = source%norm
          trget%abs_norm      = source%abs_norm
          trget%norm_sqr      = source%norm_sqr
          trget%n_tot_entries = source%n_tot_entries 

        end subroutine DS_copy_dimension

!       ----------------------------------------------------------------------
!       This subroutine removes a dimension at index dim_index from a given grid 
!       ----------------------------------------------------------------------
        subroutine DS_remove_dimension(dim_name)
        implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in) :: dim_name
!
!         Local variables
!
          integer         :: ref_dim_index, run_dim_index
!
!         Begin code
!
          ref_dim_index = DS_dim_index(ref_grid, dim_name)
          run_dim_index = DS_dim_index(run_grid, dim_name)
          call DS_remove_dimension_from_grid(ref_grid, ref_dim_index)
          call DS_remove_dimension_from_grid(run_grid, run_dim_index)
        end subroutine DS_remove_dimension
!

!       ----------------------------------------------------------------------
!       This subroutine removes a dimension at index dim_index from a given grid 
!       ----------------------------------------------------------------------
        subroutine DS_remove_dimension_from_grid(grid, dim_index)
        implicit none
!         
!         Subroutine arguments
!
          type(sampledDimension), dimension(:), allocatable,
     &      intent(inout)                          :: grid
          integer, intent(in)                      :: dim_index
!
!         Local variables
!
          type(sampledDimension), dimension(:), allocatable :: tmp
          integer i
!
!         Begin code
!

          allocate(tmp(size(grid)-1))
          do i=1,dim_index-1
            call DS_copy_dimension(grid(i), tmp(i))
          enddo
          do i=dim_index+1,size(grid)
            call DS_copy_dimension(grid(i), tmp(i-1))
          enddo
          call DS_deallocate_grid(grid)
          allocate(grid(size(tmp)))
          do i=1,size(tmp)
            call DS_copy_dimension(tmp(i), grid(i))
          enddo
          call DS_deallocate_grid(tmp)
        end subroutine DS_remove_dimension_from_grid

!       ---------------------------------------------------------------
!       This subroutine takes care of reinitializing a given dimension
!       with default values
!       ---------------------------------------------------------------
        subroutine DS_reinitialize_dimension(d_dim)
        implicit none
!         
!         Subroutine arguments
!
          type(sampledDimension), intent(inout) :: d_dim 
!
!         Local variables
!

          integer i
!
!         Begin code
!
          do i=1, size(d_dim%bins)
            call DS_reinitialize_bin(d_dim%bins(i))
          enddo
          d_dim%norm_sqr        = 0.0d0
          d_dim%abs_norm        = 0.0d0          
          d_dim%norm            = 0.0d0
          d_dim%n_tot_entries   = 0

        end subroutine DS_reinitialize_dimension

!       ---------------------------------------------------------------
!       This subroutine takes care of initializing a given dimension
!       with default values
!       ---------------------------------------------------------------
        subroutine DS_initialize_dimension(d_dim)
        implicit none
!         
!         Subroutine arguments
!
          type(sampledDimension), intent(inout) :: d_dim 
!
!         Local variables
!

          integer i
!
!         Begin code
!
          call DS_reinitialize_dimension(d_dim)
          do i=1, size(d_dim%bins)
            call DS_initialize_bin(d_dim%bins(i))
          enddo
          do i= 1, len(d_dim%dimension_name)
            d_dim%dimension_name(i:i) = ' '
          enddo
!         By default give sequential ids to the bins
          do i=1, size(d_dim%bins)
            d_dim%bins(i)%bid = i
          enddo
        end subroutine DS_initialize_dimension

!       ---------------------------------------------------------------
!       This subroutine takes care of reinitializing a given bin 
!       ---------------------------------------------------------------
        subroutine DS_initialize_bin(d_bin)
        implicit none
!         
!         Subroutine arguments
!
          type(bin), intent(inout) :: d_bin
!
!         Begin code
!
          call DS_reinitialize_bin(d_bin)
          d_bin%bid         = 0
        end subroutine DS_initialize_bin

!       ---------------------------------------------------------------
!       This subroutine takes care of initializing a given bin 
!       ---------------------------------------------------------------
        subroutine DS_reinitialize_bin(d_bin)
        implicit none
!         
!         Subroutine arguments
!
          type(bin), intent(inout) :: d_bin
!
!         Begin code
!
          d_bin%weight_sqr = 0.0d0
          d_bin%abs_weight = 0.0d0          
          d_bin%weight     = 0.0d0
          d_bin%n_entries  = 0
        end subroutine DS_reinitialize_bin

        function DS_get_dimension(grid, dim_name)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable
     &                                  :: grid
          character(len=*), intent(in)  :: dim_name
          type(sampledDimension)        :: DS_get_dimension
!
!         Begin code
!
          DS_get_dimension = grid(DS_dim_index(grid,dim_name))
        end function DS_get_dimension

!       ---------------------------------------------------------------
!       Returns the index of a bin with mBinID in the list bins
!       ---------------------------------------------------------------
        function DS_bin_index_default(bins, mBinID)
        implicit none
!         
!         Function arguments
!
          type(Bin), dimension(:), intent(in)  
     &                                  :: bins
          type(BinID)                   :: mBinID
          integer                       :: DS_bin_index_default
!
!         Begin code
!
          DS_bin_index_default = DS_bin_index_with_force(bins,mBinID,
     &                                                          .False.)
        end function DS_bin_index_default

        function DS_bin_index_with_force(bins, mBinID,force)
        implicit none
!         
!         Function arguments
!
          type(Bin), dimension(:), intent(in)  
     &                                  :: bins
          type(BinID)                   :: mBinID
          integer                       :: DS_bin_index_with_force
          logical                       :: force
!
!         Local variables
!
          integer i
!
!         Begin code
!
!         For efficiency first look at index mBinID%id
          if (size(bins).ge.mBinID%id) then
            if (bins(mBinID%id)%bid==mBinID) then
              DS_bin_index_with_force = mBinID%id
              return
            endif
          endif
          
          DS_bin_index_with_force = -1
          do i = 1, size(bins)
            if (bins(i)%bid==mBinID) then
              DS_bin_index_with_force = i
              return              
            endif
          enddo
          if (DS_bin_index_with_force.eq.-1.and.(.not.Force)) then
            write(*,*) 'DiscreteSampler:: Error in function bin_index'//
     &        "(), bin with BinID '"//trim(DS_toStr(mBinID))
     &        //"' not found."
            stop 1
          endif
        end function DS_bin_index_with_force

!       ---------------------------------------------------------------
!       Functions of the interface get_bin facilitating the access to a
!       given bin.
!       ---------------------------------------------------------------
        
        function DS_get_bin_from_binID(bins, mBinID)
        implicit none
!         
!         Function arguments
!
          type(Bin), dimension(:), intent(in)  
     &                                  :: bins
          type(BinID)                   :: mBinID
          type(Bin)                     :: DS_get_bin_from_binID
!
!         Local variables
!
          integer i
!
!         Begin code
!
          DS_get_bin_from_binID = bins(DS_bin_index(bins,mBinID))
        end function DS_get_bin_from_binID

        function DS_get_bin_from_binID_and_dimName(grid, dim_name,
     &                                                          mBinID)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable
     &                                  :: grid
          character(len=*), intent(in)  :: dim_name
          type(BinID)                   :: mBinID
          type(Bin)             :: DS_get_bin_from_binID_and_dimName
!
!         Local variables
!
          integer i
          type(SampledDimension)        :: m_dim
!
!         Begin code
!
          m_dim = DS_get_dimension(grid,dim_name)
          DS_get_bin_from_binID_and_dimName = DS_get_bin_from_binID(
     &                  m_dim%bins,mBinID)
        end function DS_get_bin_from_binID_and_dimName


!       ---------------------------------------------------------------
!       Add a new weight to a certan bin (characterized by either its 
!       binID or index)
!       ---------------------------------------------------------------
        subroutine DS_add_entry_with_BinID(dim_name, mBinID,weight)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)  :: dim_name
          type(BinID)                   :: mBinID
          real*8                        :: weight
!
!         Local variables
!
          integer dim_index, bin_index
          type(Bin)                     :: newBin
!
!         Begin code
!
          dim_index = DS_dim_index(run_grid, dim_name, .TRUE.)
          if (dim_index.eq.-1) then
              call DS_Logger('Dimension  '//dim_name//
     &        ' does not exist in the run grid. Creating it now.')
              call DS_register_dimension(dim_name,0)
              dim_index = DS_dim_index(run_grid, dim_name)
          endif

          bin_index = DS_bin_index(
     &                           run_grid(dim_index)%bins,mBinID,.TRUE.)
          if (bin_index.eq.-1) then
              call DS_Logger('Bin with binID '//trim(DS_toStr(mBinID))//
     &        ' does not exist in the run grid. Creating it now.')
              call DS_reinitialize_bin(newBin)
              newBin%bid = mBinID
              call DS_add_bin_to_bins(run_grid(dim_index)%bins,newBin)
              bin_index = DS_bin_index(run_grid(dim_index)%bins,mBinID)
          endif

!         Update global cumulative information in the grid
          run_grid(dim_index)%norm = 
     &                   run_grid(dim_index)%norm + weight
          run_grid(dim_index)%norm_sqr = 
     &                      run_grid(dim_index)%norm_sqr + weight**2
          run_grid(dim_index)%abs_norm = 
     &                    run_grid(dim_index)%abs_norm + abs(weight)
          run_grid(dim_index)%n_tot_entries =  
     &                   run_grid(dim_index)%n_tot_entries + 1
!         Update the information directly stored in the bin
          run_grid(dim_index)%bins(bin_index)%weight = 
     &           run_grid(dim_index)%bins(bin_index)%weight + weight
          run_grid(dim_index)%bins(bin_index)%weight_sqr = 
     &       run_grid(dim_index)%bins(bin_index)%weight_sqr + 
     &                                                       weight**2
          run_grid(dim_index)%bins(bin_index)%abs_weight = 
     &       run_grid(dim_index)%bins(bin_index)%abs_weight + 
     &                                                       abs(weight)
          run_grid(dim_index)%bins(bin_index)%n_entries =
     &             run_grid(dim_index)%bins(bin_index)%n_entries + 1
        end subroutine DS_add_entry_with_BinID

        subroutine DS_add_entry_with_BinIntID(dim_name, BinIntID,
     &                                                       weight)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)  :: dim_name
          integer                       :: BinIntID
          real*8                        :: weight 
!
!         Begin code
!
          call DS_add_entry_with_BinID(dim_name, DS_BinID(BinIntID),
     &                                                          weight)
        end subroutine DS_add_entry_with_BinIntID

!       ---------------------------------------------------------------
!       Prints out all informations for dimension of index d_index, or
!       name d_name.
!       ---------------------------------------------------------------
        subroutine DS_print_dim_global_info_from_void()
          integer i
          if(allocated(ref_grid).and.allocated(run_grid)) then
            do i = 1, size(ref_grid)
              call DS_print_dim_global_info_from_name(
     &                          trim(toStr(ref_grid(i)%dimension_name)))
            enddo
          else
            write(*,*) 'DiscreteSampler:: No dimension setup yet.'
          endif
        end subroutine DS_print_dim_global_info_from_void

        subroutine DS_print_dim_global_info_from_name(d_name)
        implicit none

!         Function arguments
!
          character(len=*), intent(in) :: d_name
!
!         Local variables
!
          integer n_bins, ref_dim_index, run_dim_index
!
!         Begin code
!
!         The running grid and ref grid must have the same number of
!         bins at this stage

          if(.not.(allocated(ref_grid).and.allocated(run_grid))) then
            write(*,*) 'DiscreteSampler:: No dimension setup yet.'
            return
          endif

          ref_dim_index = DS_dim_index(ref_grid,d_name,.TRUE.)
          run_dim_index = DS_dim_index(run_grid,d_name,.TRUE.)

          n_bins = size(ref_grid(DS_dim_index(ref_grid,d_name))%bins)
          write(*,*) "DiscreteSampler:: ========================"//
     &       "=========================="
          write(*,*) "DiscreteSampler:: Information for dimension '"//
     &                     d_name//"' ("//trim(toStr(n_bins))//" bins):"
          if (ref_dim_index.ne.-1) then
            write(*,*) "DiscreteSampler:: || Reference grid "
            call DS_print_dim_info(ref_grid(ref_dim_index))
          else
            write(*,*) "DiscreteSampler:: || No reference grid for "//
     &         "that dimension."
          endif
          if (run_dim_index.ne.-1) then
            write(*,*) "DiscreteSampler:: || Running grid "
            call DS_print_dim_info(run_grid(run_dim_index))
          else
            write(*,*) "DiscreteSampler:: || No running grid for "//
     &         "that dimension."
          endif
          write(*,*) "DiscreteSampler:: ========================"//
     &       "=========================="
        end subroutine DS_print_dim_global_info_from_name

!       ---------------------------------------------------------------
!       Print all informations related to a specific sampled dimension
!       in a given grid
!       ---------------------------------------------------------------
        subroutine DS_print_dim_info(d_dim)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), intent(in)  :: d_dim
!
!         Local variables
!
          integer i,j, curr_pos1, curr_pos2
          integer n_bins, bin_width
!         Adding the minimum size for the separators '|' and binID assumed
!         of being of length 2 at most, so 10*2+11 and + 20 security :)

          character(samplingBarWidth+10*2+11+20)       :: samplingBar1
          character(samplingBarWidth+10*2+11+20)       :: samplingBar2
!
!         Begin code
!
!
!         Setup the sampling bars
!
          if (d_dim%norm.eq.0.0d0) then
            samplingBar1 = "| Empty grid |"
            samplingBar2 = "| Empty grid |"
          else
            do i=1,len(samplingBar1)
              samplingBar1(i:i)=' '
              samplingBar2(i:i)=' '
            enddo
            samplingBar1(1:1) = '|'
            samplingBar2(1:1) = '|' 
            curr_pos1 = 2
            curr_pos2 = 2
            do i=1,min(10,size(d_dim%bins)) 
              samplingBar1(curr_pos1:curr_pos1+1) =
     &                             trim(DS_toStr(d_dim%bins(i)%bid))
              samplingBar2(curr_pos2:curr_pos2+1) = 
     &                             trim(DS_toStr(d_dim%bins(i)%bid))
              curr_pos1 = curr_pos1+2
              curr_pos2 = curr_pos2+2

              bin_width = int((d_dim%bins(i)%abs_weight/d_dim%abs_norm)*
     &                                                 samplingBarWidth)
              do j=1,bin_width
                samplingBar1(curr_pos1+j:curr_pos1+j) = ' '
              enddo
              curr_pos1 = curr_pos1+bin_width+1
              samplingBar1(curr_pos1:curr_pos1) = '|'
              curr_pos1 = curr_pos1+1

              bin_width = int((float(d_dim%bins(i)%n_entries)/
     &                            d_dim%n_tot_entries)*samplingBarWidth)
              do j=1,bin_width
                samplingBar2(curr_pos2+j:curr_pos2+j) = ' '
              enddo
              curr_pos2 = curr_pos2+bin_width+1
              samplingBar2(curr_pos2:curr_pos2) = '|'
              curr_pos2 = curr_pos2+1
            enddo
          endif
!
!         Write out info
!
          n_bins = size(d_dim%bins)
          
          write(*,*) "DiscreteSampler::   -> Total number of "//
     &         "entries : "//trim(toStr(d_dim%n_tot_entries))
          if (n_bins.gt.10) then
            write(*,*) "DiscreteSampler::   -> Sampled as"//
     &                                      " (first 10 bins):"
          else
            write(*,*) "DiscreteSampler::   -> Sampled as:"
          endif
          write(*,*) "DiscreteSampler::    "//trim(samplingBar2)
          write(*,*) "DiscreteSampler::   -> (norm_sqr , "//
     &      "abs_norm , norm     ) :"
          write(*,*) "DiscreteSampler::      ("//
     &      trim(toStr(d_dim%norm_sqr,'Ew.3'))//", "//
     &      trim(toStr(d_dim%abs_norm,'Ew.3'))//", "//
     &      trim(toStr(d_dim%norm,'Ew.3'))//")"
          if (n_bins.gt.10) then
            write(*,*) "DiscreteSampler::   -> Sampled as"//
     &                                      " (first 10 bins):"
          else
            write(*,*) "DiscreteSampler::   -> Sampled as:"
          endif
          write(*,*) "DiscreteSampler::    "//trim(samplingBar1)

        end subroutine DS_print_dim_info

!       ---------------------------------------------------------------
!         Functions to add a bin with different binID specifier
!       ---------------------------------------------------------------      
        subroutine DS_add_bin_with_IntegerID(dim_name,intID)
          implicit none
!         
!         Subroutine arguments
!
          integer, intent(in)      :: intID
          character(len=*)         :: dim_name
!
!         Begin code
!
          call DS_add_bin_with_binID(dim_name,DS_binID(intID))
        end subroutine DS_add_bin_with_IntegerID

        subroutine DS_add_bin_with_void(dim_name)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*)         :: dim_name
!
!         Local variables
!
          integer                  :: ref_size, run_size
!
!         Begin code
!
          ref_size=size(ref_grid(DS_dim_index(ref_grid,dim_name))%bins)
          run_size=size(run_grid(DS_dim_index(run_grid,dim_name))%bins)
          call DS_add_bin_with_binID(dim_name,DS_binID(
     &                                max(ref_size, run_size)+1))
        end subroutine DS_add_bin_with_void

        subroutine DS_add_bin_with_binID(dim_name,mBinID)
          implicit none
!         
!         Subroutine arguments
!
          type(binID), intent(in)  :: mBinID
          character(len=*)         :: dim_name
!
!         Local variables
!
          type(Bin)                :: new_bin
!
!         Begin code
!
          call DS_reinitialize_bin(new_bin)
          new_bin%bid = mBinID
          call DS_add_bin_to_bins(ref_grid(DS_dim_index(ref_grid, 
     &                                          dim_name))%bins,new_bin)
          call DS_add_bin_to_bins(run_grid(DS_dim_index(run_grid, 
     &                                          dim_name))%bins,new_bin)
        end subroutine DS_add_bin_with_binID

        subroutine DS_add_bin_to_bins(bins,new_bin)
          implicit none
!         
!         Subroutine arguments
!
          type(Bin), dimension(:), allocatable, intent(inout)  
     &                             :: bins
          type(Bin)                :: new_bin
!
!         Local variables
!
          type(Bin), dimension(:), allocatable :: tmp
          integer                              :: i, bin_index
!
!         Begin code
!
          bin_index = DS_bin_index(bins,new_bin%bid,.True.)
          if (bin_index.ne.-1) then
             write(*,*)"DiscreteSampler:: Error, the bin with binID '"//
     &         trim(DS_toStr(new_bin%bid))//"' cannot be added "//
     &         "be added because it already exists."
               stop 1
          endif


          allocate(tmp(size(bins)+1))
          do i=1,size(bins)
            call DS_copy_bin(bins(i),tmp(i))
          enddo
          tmp(size(bins)+1) = new_bin
          deallocate(bins)
          allocate(bins(size(tmp)))
          do i=1,size(bins)
            call DS_copy_bin(tmp(i),bins(i))          
          enddo
          deallocate(tmp)
        end subroutine DS_add_bin_to_bins

        subroutine DS_copy_bin(source, trget)
            implicit none
            type(Bin), intent(out) :: trget
            type(Bin), intent(in)  :: source
            trget%weight     = source%weight
            trget%weight_sqr = source%weight_sqr
            trget%abs_weight = source%abs_weight
            trget%n_entries  = source%n_entries
            trget%bid        = DS_binID(source%bid%id)
        end subroutine DS_copy_bin

!       ---------------------------------------------------------------
!         Functions to remove a bin from a dimension
!       ---------------------------------------------------------------
        subroutine DS_remove_bin_withIndex(dim_name, binIndex)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)   :: dim_name
          integer, intent(in)            :: binIndex
!
!         Begin code
!

          call DS_remove_bin_from_grid(run_grid(
     &                       DS_dim_index(run_grid, dim_name)),binIndex)
        end subroutine DS_remove_bin_withIndex

        subroutine DS_remove_bin_withBinID(dim_name, mbinID)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)   :: dim_name
          type(binID), intent(in)        :: mbinID
!
!         Local variables
!
          integer                        :: ref_dim_index,run_dim_index
          integer                        :: ref_bin_index,run_bin_index
!
!         Begin code
!
          ref_dim_index = DS_dim_index(ref_grid, dim_name)
          ref_bin_index = DS_bin_index(ref_grid(ref_dim_index)%bins,
     &                                                          mbinID)
          call DS_remove_bin_from_grid(ref_grid(ref_dim_index),
     &                                                   ref_bin_index)
          run_dim_index = DS_dim_index(run_grid, dim_name)
          run_bin_index = DS_bin_index(run_grid(run_dim_index)%bins,
     &                                                          mbinID)
          call DS_remove_bin_from_grid(run_grid(run_dim_index),
     &                                                   run_bin_index)
        end subroutine DS_remove_bin_withBinID

        subroutine DS_remove_bin_withIntegerID(dim_name, mBinIntID)
          implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)   :: dim_name
          integer, intent(in)            :: mBinIntID       
!
!         Begin code
!
          call DS_remove_bin_withBinID(dim_name,DS_binID(mBinIntID))
        end subroutine DS_remove_bin_withIntegerID

        subroutine DS_remove_bin_from_grid(grid, bin_index)
          implicit none
!         
!         Subroutine arguments
!
          type(SampledDimension), intent(inout)  :: grid
          integer, intent(in)                    :: bin_index
!
!         Local variables
!
          type(Bin), dimension(:), allocatable :: tmp
          integer                              :: i
!
!         Begin code
!

!         Update the norm, norm_sqr and the number of entries in
!         the corresponding dimension
          grid%norm = grid%norm - grid%bins(bin_index)%weight
          grid%norm_sqr = grid%norm_sqr - 
     &                                   grid%bins(bin_index)%weight_sqr
          grid%abs_norm = grid%abs_norm -
     &                                   grid%bins(bin_index)%abs_weight
          grid%n_tot_entries = grid%n_tot_entries
     &                                  - grid%bins(bin_index)%n_entries

          allocate(tmp(size(grid%bins)-1))
          do i=1,bin_index-1
            tmp(i) = grid%bins(i)
          enddo
          do i=bin_index+1,size(grid%bins)
            tmp(i-1) = grid%bins(i)          
          enddo
          deallocate(grid%bins)
          allocate(grid%bins(size(tmp)))
          do i=1,size(tmp)
            grid%bins(i)=tmp(i)
          enddo
          deallocate(tmp)
        end subroutine DS_remove_bin_from_grid


!       ---------------------------------------------------------------
!       Function to update the reference grid with the running one
!       ---------------------------------------------------------------
        subroutine DS_update_all_grids(filterZeros)
        implicit none
!         
!         Subroutine arguments
!
          logical, optional :: filterZeros
!         
!         Local variables
!
          integer           :: i
          logical           :: do_filterZeros          
!
!         Begin code
!
          if(present(filterZeros)) then
            do_filterZeros = filterZeros
          else
            do_filterZeros = .False.
          endif
          do i=1, size(run_grid)
            call DS_update_grid_with_dim_index(i,do_filterZeros)
          enddo
        end subroutine DS_update_all_grids

        subroutine DS_update_grid_with_dim_name(dim_name, filterZeros)
        implicit none
!         
!         Subroutine arguments
!
          character(len=*)                 :: dim_name
          logical, optional                :: filterZeros          
!         
!         Local variables
!
          integer           :: i
          logical           :: do_filterZeros
!
!         Begin code
!
          if(present(filterZeros)) then
            do_filterZeros = filterZeros
          else
            do_filterZeros = .False.
          endif
          call DS_update_grid_with_dim_index(
     &                   DS_dim_index(run_grid,dim_name),do_filterZeros)

        end subroutine DS_update_grid_with_dim_name

        subroutine DS_update_grid_with_dim_index(d_index,filterOutZeros)
        implicit none
!         
!         Subroutine arguments
!
          integer                               :: d_index
          logical                               :: filterOutZeros
!         
!         Local variables
!
          integer                               :: i, ref_d_index 
          integer                               :: ref_bin_index
          integer                               :: j, shift
          character, dimension(:), allocatable  :: dim_name
          type(BinID)                           :: mBinID
          type(Bin)                             :: new_bin
!
!         Begin code
!
          allocate(dim_name(size(run_grid(d_index)%dimension_name)))
          dim_name = run_grid(d_index)%dimension_name
          call DS_Logger("Updating dimension '"//
     &                                      trim(toStr(dim_name))//"'.")

!         Start by making sure that the dimension exists in the
!         reference grid. If not, then create it.
          if (DS_dim_index(ref_grid,
     &         run_grid(d_index)%dimension_name,.True.).eq.-1) then
              call DS_Logger('Reference grid does not have dimension '//
     &                         trim(toStr(dim_name))//'. Adding it now')
              call DS_add_dimension_to_grid(ref_grid, 
     &                                        trim(toStr(dim_name)) , 0)
          endif

          ref_d_index = DS_dim_index(ref_grid, dim_name)

          do i=1,size(run_grid(d_index)%bins)
            mBinID = run_grid(d_index)%bins(i)%bid
            ref_bin_index = DS_bin_index(
     &                        ref_grid(ref_d_index)%bins,mBinID,.True.) 
            if (ref_bin_index.eq.-1) then
              call DS_Logger('Bin with binID '//trim(DS_toStr(mBinID))//
     &              ' is missing in the reference grid. Adding it now.')
              call DS_reinitialize_bin(new_bin)
              new_bin%bid = mBinID
              call DS_add_bin_to_bins(ref_grid(ref_d_index)%bins,
     &                                                          new_bin)
              ref_bin_index = DS_bin_index(
     &                                ref_grid(ref_d_index)%bins,mBinID)
            endif

            new_bin = ref_grid(ref_d_index)%bins(ref_bin_index) +
     &                                         run_grid(d_index)%bins(i)

            ref_grid(ref_d_index)%bins(ref_bin_index) = new_bin
            ref_grid(ref_d_index)%norm = 
     &                       ref_grid(ref_d_index)%norm + new_bin%weight
            ref_grid(ref_d_index)%norm_sqr = 
     &               ref_grid(ref_d_index)%norm_sqr + new_bin%weight_sqr
            ref_grid(ref_d_index)%abs_norm = 
     &               ref_grid(ref_d_index)%abs_norm + new_bin%abs_weight
            ref_grid(ref_d_index)%n_tot_entries = 
     &           ref_grid(ref_d_index)%n_tot_entries + new_bin%n_entries
          enddo

!         Now filter all bins in ref_grid that have 0.0 weight and
!         remove them! They will not be probed anyway.
          if (filterOutZeros) then
            shift = 0
            do j=1,size(ref_grid(ref_d_index)%bins)
              i = j - shift
              if ((ref_grid(ref_d_index)%bins(i)%weight.eq.0.0d0).and.
     &        (ref_grid(ref_d_index)%bins(i)%abs_weight.eq.0.0d0).and.
     &        (ref_grid(ref_d_index)%bins(i)%weight_sqr.eq.0.0d0)) then
                call DS_Logger('Bin with binID '//
     &            trim(DS_toStr(ref_grid(ref_d_index)%bins(i)%bid))//
     &            ' is zero and will be filtered out. Removing it now.')
                call DS_remove_bin_from_grid(ref_grid(ref_d_index),i)
                shift = shift + 1
              endif
            enddo
          endif

!         Clear the running grid now
          call DS_reinitialize_dimension(run_grid(d_index))

          deallocate(dim_name)

        end subroutine DS_update_grid_with_dim_index


        function DS_combine_two_bins(BinA, BinB) result(CombinedBin)
        implicit none
!         
!         Function arguments
!
          integer               :: d_index
          Type(Bin), intent(in) :: BinA, BinB
          Type(Bin)             :: CombinedBin
!         
!         Local variables
!
          call DS_reinitialize_bin(CombinedBin)
          if(.not.(BinA%bid==BinB%bid)) then
            write(*,*) 'DiscreteSampler:: Error in function '//
     &        'DS_combine_two_bins, cannot combine two bins '//
     &        ' with different bin IDs : '//trim(DS_toStr(BinA%bid))//
     &        ', '//trim(DS_toStr(BinB%bid))
            stop 1
          endif
          CombinedBin%bid = BinA%bid
          CombinedBin%weight = BinA%weight + BinB%weight
          CombinedBin%abs_weight = BinA%abs_weight + BinB%abs_weight
          CombinedBin%weight_sqr = BinA%weight_sqr + BinB%weight_sqr
          CombinedBin%n_entries = BinA%n_entries + BinB%n_entries
        end function DS_combine_two_bins

!       ================================================
!       Grid I/O functions
!       ================================================

!       ---------------------------------------------------------------
!       This function writes the ref_grid to a file specified by its 
!       filename.
!       ---------------------------------------------------------------
        subroutine DS_write_grid_with_filename(filename, dim_name)
        implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)           :: filename
          character(len=*), intent(in), optional :: dim_name
!         
!         Local variables
!
          logical fileExist
!
!         Begin code
!
          inquire(file=filename, exist=fileExist)
          if (fileExist) then
            call DS_Logger('DiscreteSampler:: The file '
     &        //filename//' already exists, so beware that '//
     &                  ' the grid information will be appended to it.')
          endif
          open(123, file=filename, err=11, access='append',
     &                                                   action='write')
          goto 12
11        continue
          write(*,*) 'DiscreteSampler :: Error, file '//filename//
     &                               ' could not be opened for writing.'
          stop 1
12        continue
          if (present(dim_name)) then
            call DS_write_grid_with_streamID(123, dim_name)
          else
            call DS_write_grid_with_streamID(123)
          endif
          close(123)
        end subroutine DS_write_grid_with_filename

!       ---------------------------------------------------------------
!       This function writes the ref_grid or all grids to a file
!       specified by its stream ID.
!       ---------------------------------------------------------------
        subroutine DS_write_grid_with_streamID(streamID, dim_name)
        implicit none
!         
!         Subroutine arguments
!
          integer, intent(in)                    :: streamID
          character(len=*), intent(in), optional :: dim_name
!         
!         Local variables
!
          type(SampledDimension)                 :: grid
          integer                                :: i
!
!         Begin code
!
          if (present(dim_name)) then
            grid = ref_grid(DS_dim_index(ref_grid, dim_name))
            call DS_write_grid_from_grid(grid, streamID)
          else
            do i=1,size(ref_grid)
              grid = ref_grid(i)
              call DS_write_grid_from_grid(grid, streamID)
            enddo
          endif
        end subroutine DS_write_grid_with_streamID

!       ---------------------------------------------------------------
!       This function writes a given grid to a file.
!       ---------------------------------------------------------------
        subroutine DS_write_grid_from_grid(grid, streamID)
        implicit none
!         
!         Subroutine arguments
!
          integer, intent(in)                    :: streamID
          type(SampledDimension), intent(in)     :: grid
!         
!         Local variables
!
          integer                                :: i
!
!         Begin code
!

          write(streamID,*) ' <DiscreteSampler_grid>'
          write(streamID,*) ' '//trim(toStr(grid%dimension_name))
          write(streamID,*) '# binID   n_entries weight   weight_sqr'//
     &      '   abs_weight'
          do i=1,size(grid%bins)
            write(streamID,*) 
     &                 '   '//trim(DS_toStr(grid%bins(i)%bid))//
     &                 '   '//trim(toStr(grid%bins(i)%n_entries))//
     &                 '   '//trim(toStr(grid%bins(i)%weight))//
     &                 '   '//trim(toStr(grid%bins(i)%weight_sqr))//
     &                 '   '//trim(toStr(grid%bins(i)%abs_weight))
          enddo
          write(streamID,*) ' </DiscreteSampler_grid>'

        end subroutine DS_write_grid_from_grid

!       ---------------------------------------------------------------
!       This function loads the grid specified in a file specified by its
!       stream ID into the run_grid.
!       ---------------------------------------------------------------
        subroutine DS_load_grid_with_filename(filename, dim_name)
        implicit none
!         
!         Subroutine arguments
!
          character(len=*), intent(in)           :: filename
          character(len=*), intent(in), optional :: dim_name
!         
!         Local variables
!
          logical fileExist        
!
!         Begin code
!
          inquire(file=filename, exist=fileExist)
          if (.not.fileExist) then
            write(*,*) 'DiscreteSampler:: Error, the file '//filename//
     &                                           ' could not be found.'
            stop 1
          endif
          open(124, file=filename, err=13, action='read')
          goto 14
13        continue
          write(*,*) 'DiscreteSampler :: Error, file '//filename//
     &                                ' exists but could not be read.'
14        continue
          if (present(dim_name)) then
            call DS_load_grid_with_streamID(124, dim_name)
          else
            call DS_load_grid_with_streamID(124)
          endif
          close(124)
        end subroutine DS_load_grid_with_filename

!       ---------------------------------------------------------------
!       This function loads the grid specified in a file specified by its 
!       stream ID into the run_grid.
!       ---------------------------------------------------------------
        subroutine DS_load_grid_with_streamID(streamID, dim_name)
        implicit none
!         
!         Subroutine arguments
!
          integer, intent(in)                    :: streamID
          character(len=*), intent(in), optional :: dim_name
!         
!         Local variables
!
          integer                                :: i
          character(512)                         :: buff
          character(2)                           :: TwoBuff
          character(3)                           :: ThreeBuff
          logical                                :: startedGrid
          real*8                       :: weight, abs_weight, weight_sqr
          integer                      :: n_entries, bid
          type(Bin)                    :: new_bin
          integer                      :: char_size

!
!         Begin code
!
!         First reinitialize the running_grid
          call DS_deallocate_grid(run_grid)
!         Now start reading the file
          startedGrid = .False.
          do
            read(streamID, "(A)", size=char_size, end=999, advance='no')
     &        TwoBuff
            if (char_size.le.1) then
              cycle
            endif
            if (TwoBuff(1:1).eq.'#'.or.TwoBuff(2:2).eq.'#') then
!             Advance the stream
              read(streamID,*,end=990) buff
              cycle
            endif
            if (startedGrid) then
              read(streamID, "(A)", size=char_size,
     &                                    end=999, advance='no') TwoBuff
              if (TwoBuff(1:2).eq.'</') then
!             Advance the stream
                read(streamID,*,end=990) buff
                startedGrid = .False.
                cycle
              endif
              read(streamID,*,end=990) bid, n_entries, weight, 
     &                                            weight_sqr, abs_weight
              new_bin%bid           = bid
              new_bin%n_entries     = n_entries
              new_bin%weight        = weight
              new_bin%weight_sqr    = weight_sqr
              new_bin%abs_weight    = abs_weight
              call DS_add_bin_to_bins(run_grid(size(run_grid))%bins,
     &                                                          new_bin)
            else
!             Advance the stream
              read(streamID,*,end=990) buff
              if (buff(1:22).eq.'<DiscreteSampler_grid>') then
                startedGrid = .True.
                read(streamID,*,end=990) buff
                call DS_register_dimension(trim(buff),0,.False.)
              endif
            endif
          enddo
          goto 999
990       continue
          write(*,*) 'DiscreteSampler:: Error, when loading grids'//
     &      ' from file.'
          stop 1
999       continue

!         Now update the running grid into the reference one
          call DS_update_grid()
        end subroutine DS_load_grid_with_streamID


!       ---------------------------------------------------------------
!       Synchronizes the cumulative information in a given grid from
!       its bins.
!       --------------------------------------------------------------- 
        subroutine DS_synchronize_grid_with_bins(grid)
        implicit none
!
!         Subroutine argument
!
          type(sampledDimension)                 :: grid
!         
!         Local variables
!
          real*8           :: norm, abs_norm, norm_sqr
          integer          :: i, n_tot_entries
!
!         Begin Code
!
          norm              = 0.0d0
          abs_norm          = 0.0d0
          norm_sqr          = 0.0d0
          n_tot_entries     = 0
          do i=1,size(grid%bins)
            n_tot_entries   = n_tot_entries  + grid%bins(i)%n_entries
            norm_sqr        = norm_sqr       + grid%bins(i)%weight_sqr
            abs_norm        = abs_norm       + grid%bins(i)%abs_weight
            norm            = norm           + grid%bins(i)%weight
          enddo
          grid%n_tot_entries = n_tot_entries
          grid%norm_sqr      = norm_sqr
          grid%abs_norm      = abs_norm
          grid%norm          = norm
        end subroutine DS_synchronize_grid_with_bins

!       ================================================
!       Functions and subroutine handling derived types
!       ================================================

!       ---------------------------------------------------------------
!       Specify how bin idea should be compared
!       ---------------------------------------------------------------
        function equal_binID(binID1,binID2)
        implicit none
!         
!         Function arguments
!
          type(binID), intent(in)  :: binID1, binID2
          logical                  :: equal_binID
!
!         Begin code
!
          if(binID1%id.ne.binID2%id) then
            equal_binID = .False.
            return
          endif
          equal_binID = .True.
          return
        end function equal_binID

!       ---------------------------------------------------------------
!       BinIDs constructors
!       ---------------------------------------------------------------
        pure elemental subroutine binID_from_binID(binID1,binID2)
        implicit none
!         
!         Function arguments
!
          type(binID), intent(out)  :: binID1
          type(binID), intent(in)  :: binID2
!
!         Begin code
!
          binID1%id = binID2%id
        end subroutine binID_from_binID

        pure elemental subroutine binID_from_integer(binID1,binIDInt)
        implicit none
!         
!         Function arguments
!
          type(binID), intent(out)  :: binID1
          integer,     intent(in)   :: binIDInt
!
!         Begin code
!
          binID1%id = binIDInt
        end subroutine binID_from_integer

!       Provide a constructor-like way of creating a binID
        function DS_binID(binIDInt)
        implicit none
!         
!         Function arguments
!
          type(binID)              :: DS_binID
          integer,     intent(in)  :: binIDInt
!
!         Begin code
!
          DS_binID = binIDInt
        end function DS_binID
!       ---------------------------------------------------------------
!       String representation of a binID
!       ---------------------------------------------------------------
        function DS_toStr(mBinID)
        implicit none
!         
!         Function arguments
!
          type(binID), intent(in)  :: mBinID
          character(100)           :: DS_toStr
!
!         Begin code
!
          DS_toStr = trim(toStr(mBinID%id))
        end function DS_toStr


!       ================================================
!        Access routines emulating a dictionary
!       ================================================

!       ---------------------------------------------------------------
!       Returns the index of the discrete dimension with name dim_name
!       ---------------------------------------------------------------
        function DS_dim_index_default(grid, dim_name)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable
     &                                  :: grid
          character(len=*), intent(in)  :: dim_name
          integer                       :: DS_dim_index_default
!
!         Begin code
!  
          DS_dim_index_default =
     &               DS_dim_index_with_force(grid, dim_name, .False.)
        end function DS_dim_index_default

        function DS_dim_index_with_force(grid, dim_name, force)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable
     &                                  :: grid
          character(len=*), intent(in)  :: dim_name
          integer                       :: DS_dim_index_with_force
          logical                       :: force
!
!         Local variables
!

          integer i,j
!
!         Begin code
!
          DS_dim_index_with_force = -1
          if (.not.allocated(grid)) then
            return
          endif
          do i = 1, size(grid)
            if (len(dim_name).ne.size(grid(i)%dimension_name)) cycle
            do j =1, len(dim_name)
              if(grid(i)%dimension_name(j).ne.dim_name(j:j)) then
                goto 1
              endif
            enddo
            DS_dim_index_with_force = i
            return
1           continue
          enddo
          if (DS_dim_index_with_force.eq.-1.and.(.not.force)) then
            write(*,*) 'DiscreteSampler:: Error in function dim_index'//
     &        "(), dimension name '"//dim_name//"' not found."
            stop 1
          endif
        end function DS_dim_index_with_force

        function DS_dim_index_default_with_chararray(grid, dim_name)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable 
     &                                         :: grid
          character, dimension(:), intent(in)  :: dim_name
          integer                 :: DS_dim_index_default_with_chararray
!
!         Begin code
!  
          DS_dim_index_default_with_chararray = 
     &                DS_dim_index_with_force_with_chararray(
     &                                          grid, dim_name, .False.)
        end function DS_dim_index_default_with_chararray

        function DS_dim_index_with_force_with_chararray(
     &                                            grid, dim_name, force)
        implicit none
!         
!         Function arguments
!
          type(sampledDimension), dimension(:), intent(in), allocatable
     &                                        :: grid
          character, dimension(:), intent(in) :: dim_name
          integer              :: DS_dim_index_with_force_with_chararray
          logical                             :: force
!
!         Local variables
!

          integer i,j
!
!         Begin code
!
          DS_dim_index_with_force_with_chararray = -1
          if (.not.allocated(grid)) then
            return
          endif
          do i = 1, size(grid)
            if (size(dim_name).ne.size(grid(i)%dimension_name)) cycle
            do j =1, size(dim_name)
              if(grid(i)%dimension_name(j).ne.dim_name(j)) then
                goto 1
              endif
            enddo
            DS_dim_index_with_force_with_chararray = i
            return
1           continue
          enddo
          if (DS_dim_index_with_force_with_chararray.eq.-1.and.
     &                                                (.not.force)) then
            write(*,*) 'DiscreteSampler:: Error in function dim_index'//
     &        "(), dimension name '"//dim_name//"' not found."
            stop 1
          endif
        end function DS_dim_index_with_force_with_chararray

!       End module
        end module DiscreteSampler
