      program tester

          use DiscreteSampler

          type(binID) myBinID
          real*8 mValue
          integer i

          call DS_register_dimension('Dimension1',10)
          mValue=0.0d0
          do i=1,10
            mValue=mValue+1.0d0
            call DS_add_entry('Dimension1',i,mValue)
          enddo
          call DS_update_grid()
          call DS_write_grid('grids.dsg')
          call DS_clear()
          call DS_register_dimension('Dimension2',10)
          mValue=0.0d0
          do i=1,10
            mValue=mValue+1.0d0
            call DS_add_entry('Dimension2',i,mValue)
          enddo
          call DS_update_grid()
          call DS_write_grid('grids.dsg')
          call DS_clear()
          write(*,*) 'Before grid load'
          call DS_print_global_info()          
          call DS_load_grid('grids.dsg')
          write(*,*) 'after grid load'
          call DS_print_global_info()          

      end program tester
