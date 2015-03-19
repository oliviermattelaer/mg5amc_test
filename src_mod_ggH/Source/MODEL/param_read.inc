ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision Hmass
      common /sushiMH/ Hmass

      INCLUDE '../param_card.inc'

      call init_sushi()
      MDL_MX0 = Hmass
