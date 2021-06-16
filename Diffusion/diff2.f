        program diffusion
        real rho(500,500)
        print*,"Enter the number of timesteps"
        read(5,*) nts ! number of timesteps
        open(7,file="rho")
  111   format(500(2x,e12.6))
        open(8,file="entr")
  112   format(2(2x,e12.5))
        rho=0.0
        do i=230,270 ! initial config
          do j=230,270
            rho(i,j)=1.0
          enddo
        enddo ! end initial config
        do it=1,nts ! calculation start, loop over timesteps
          s=0.0
          do i=2,499
            do j=2,499
              rho(i,j)=rho(i,j)+0.3*(rho(i+1,j)+rho(i-1,j)+rho(i,j+1)
     *        +rho(i,j-1)-4.0*rho(i,j)) ! continuation by symbol at pos 6
              rh=rho(i,j)
              if(rh>0.0) then
                s=s-rh*log(rh)
              endif
            enddo ! i
          enddo ! j
          sit=it
          write(8,112) sit,s
        enddo ! timesteps
        do i=1,500 ! print results
          write(7,111) (rho(i,j),j=1,500)
        enddo
        end







