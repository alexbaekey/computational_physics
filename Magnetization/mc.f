        program montecarlo
        real sm(10,10),esm(10,10) ! sm is spin, esm is energy
        ej=1.0 ! coupling constant
        dt=0.03 ! termperature step
        open(7,file="mag_mc")
  111   format(3(2x,e12.5))
        print*,"Enter the temperature index"
        read(5,*) itemp
        smtot=0.0
        do i=1,10
          do j=1,10
            sm(i,j)=1.0
            r1=rand(0)
            if(r1<0.5) then
              sm(i,j)=-1.0
            endif
            smtot=smtot+sm(i,j)
          enddo
        enddo
        print*, smtot
        temp=itemp*dt
        do its=1,1000 ! time steps 
        smtot=0.0
        etot=0.0 
        do i=2,9
          do j=2,9
            sm4=sm(i,j-1)+sm(i,j+1)+sm(i-1,j)+sm(i+1,j)
            esm(i,j)=-ej*sm(i,j)*sm4
            call flip(sm(i,j),esm(i,j),temp)
            etot=etot+esm(i,j)
            smtot=smtot+sm(i,j)
          enddo
        enddo
        do j=2,9 ! vertical edges
          sm1j=sm(1,j-1)+sm(1,j+1)+sm(2,j)+sm(10,j) ! periodic boundary
          esm(1,j)=-ej*sm(1,j)*sm1j
          call flip(sm(1,j),esm(1,j),temp)
          etot=etot+esm(1,j)
          smtot=smtot+sm(1,j)
          sm10j=sm(10,j-1)+sm(10,j+1)+sm(9,j)+sm(1,j) ! periodic boundary
          esm(10,j)=-ej*sm(10,j)*sm10j
          call flip(sm(10,j),esm(10,j),temp)
          etot=etot+esm(10,j)
          smtot=smtot+sm(10,j)
        enddo ! vertical edges
        do i=2,9 ! horizontal edges
          sm1i=sm(i-1,1)+sm(i+1,1)+sm(i,2)+sm(i,10) ! periodic boundary
          esm(i,1)=-ej*sm(i,1)*sm1i
          call flip(sm(i,1),esm(i,1),temp)
          etot=etot+esm(i,1)
          smtot=smtot+sm(i,1)
          sm10i=sm(i-1,10)+sm(i+1,10)+sm(i,9)+sm(i,1) ! periodic boundary
          esm(i,10)=-ej*sm(i,10)*sm10i
          call flip(sm(i,10),esm(i,10),temp)
          etot=etot+esm(i,10)
          smtot=smtot+sm(i,10)
        enddo ! horizontal edges 
        !corners
        !bottom left 1,1
        sm11=sm(1,2)+sm(2,1)+sm(1,10)+sm(10,1)
        esm(1,1)=-ej*sm(1,1)*sm11
        call flip(sm(1,1),esm(1,1),temp) 
        smtot=smtot+sm(1,1)
        etot=etot+esm(1,1)
        !bottom right 10,1
        sm101=sm(10,2)+sm(9,1)+sm(1,1)+sm(10,10)
        esm(10,1)=-ej*sm(10,1)*sm101
        call flip(sm(10,1),esm(10,1),temp) 
        smtot=smtot+sm(10,1)
        etot=etot+esm(10,1)     
        !top left 1,10
        sm110=sm(1,9)+sm(2,10)+sm(1,1)+sm(10,10)
        esm(1,10)=-ej*sm(1,10)*sm110
        call flip(sm(1,10),esm(1,10),temp) 
        smtot=smtot+sm(1,10)
        etot=etot+esm(1,10)
        !top right 10,10
        sm1010=sm(10,9)+sm(9,10)+sm(1,10)+sm(10,1)
        esm(10,10)=-ej*sm(10,10)*sm1010
        call flip(sm(10,10),esm(10,10),temp) ! decides flip or not
        smtot=smtot+sm(10,10)
        etot=etot+esm(10,10)
        ! algorithm to flip or not
        ! print*, "Etot, stot",etot,smtot
        ts=its
        write(7,111) ts,smtot,etot ! time. spin, energy
        enddo ! timesteps
        end

        subroutine flip(smf,esmf,temp) ! spin, energy, temp for flip
        esm1=-esmf ! flip
        if(esm1<esmf) then
          esmf=esm1
          smf=-smf
        else
          bf=exp(-2.0*esm1/temp) ! calculate boltzman factor
          r1=rand(0)
          if(r1<bf) then
            smf=-smf
            esmf=esm1
          endif
        endif
        return 
        end






