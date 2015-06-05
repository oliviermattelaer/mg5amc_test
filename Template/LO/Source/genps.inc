c*************************************************************************
c       Parameters used by genps and dsample, you must recompile
c	dsample if you change anything below
c*************************************************************************
      include 'maxparticles.inc'
      integer    max_branch
      parameter (max_branch=max_particles-1)
      integer    ng     , maxdim     , maxinvar     , maxconfigs
      parameter (ng = 96, maxdim = 3*(max_particles-2)-1, maxinvar= 4*max_particles, maxconfigs=10)
c      parameter (ng = 43, maxdim = 25, maxinvar= 57 , maxconfigs=20)

      double precision  xgmin,     xgmax
      parameter        (xgmin=-1d0, xgmax=1d0)

      integer    maxevents         !Requires about 1K/event
      parameter (maxevents=200000)  !Maximum # events to write to disk

c*************************************************************************
c	Parameters used for parallel running
c*************************************************************************
      integer  maxplace    ,maxpoints    ,maxans
      parameter (maxplace=9,maxpoints=10,maxans=5)
c      parameter (maxplace=199,maxpoints=100,maxans=50)
      integer    maxprb
      parameter (maxprb = maxconfigs*maxplace*maxpoints)
      integer    maxfprb
      parameter (maxfprb = maxinvar*maxplace*maxpoints)
c*************************************************************************
c   NB_subchannel (number of multi-channel handle at fortran level)
c************************************************************************* 
      integer nb_channel
      parameter (nb_channel=1)
c     store which of the change of variable/integrand to use
      integer ch_picked
      common/secondary_multichannel/ch_picked

     
c*************************************************************************
c	Parameters for helicity and color flow sums in matrixN.f
c*************************************************************************
      REAL*8 LIMHEL
      PARAMETER(LIMHEL=1e-6)
      INTEGER MAXTRIES,NUMPERMS
      PARAMETER(MAXTRIES=100,NUMPERMS=10)
      INTEGER NLCFACT
      PARAMETER(NLCFACT=10)
C     To pass the helicity/ColorFlow configuration chosen by the DiscreteSampler to 
C     matrix<i>.f
      double precision hel_jacobian
      INTEGER HEL_PICKED
      COMMON/HEL_PICKED/HEL_PICKED,hel_jacobian
      double precision cf_jacobian
      INTEGER CF_PICKED
      COMMON/CF_PICKED/CF_PICKED,cf_jacobian
C    
C     PARAMETER TO DISTINGUISH Color ordered computation from the rest
      LOGICAL COLOR_ORDERED
      INTEGER MC_symfact
      COMMON/to_Discrete_sampler_mode/COLOR_ORDERED,MC_symfact
    