        program dla2 ! now includes dimenstionality analysis
        integer iclx(1000), icly(1000)
        open(7,file="fractal")
  111   format(2(2x,i4)) ! 2 clolumn, x and y, 2 spaces, 4 digits
        open(8,file="dim")
  112   format(2(2x,e12.5)) 
        print*, "Enter the final size of the cluster"
        read(5,*) nacf
        nac=1 ! current number of atoms belonging to the cluster
        iclx(1)=0
        icly(1)=0
   2    continue
        rad=0.0
        do i=1,nac
          rx=abs(iclx(i))
          if(rx>rad) then
            rad=rx
          endif
          ry=abs(icly(i))
          if(ry>rad) then
            rad=ry
          endif
        enddo
        rad=rad+5.0
        do iw=1,50000 ! number of walkers
          r1=rand(0) ! modulo of x
          atx=r1*rad
          iatx=int(atx)
          r2=rand(0)
          if(r2<0.5) then
            iatx=-iatx
          endif
          iaty=int(sqrt(rad**2-atx**2))
          r3=rand(0)
          if(r3<0.5) then
            iaty=-iaty
          endif
          do it=1,10000 !timesteps
            r2=rand(0)
            if(r2<0.5) then ! along x
              if(r2<0.25) then 
                iatx=iatx+1
              else
                iatx=iatx-1
              endif
            else ! along y
              if(r2<0.75) then
                iaty=iaty+1
              else
                iaty=iaty-1
              endif
            endif
            do i=1,nac ! is it a neighbor?
              difx=iatx-iclx(i)
              dify=iaty-icly(i)
              dif=sqrt(difx**2+dify**2)
              if(dif<1.001) then
                goto 1
              endif
            enddo
          enddo ! timesteps
        enddo ! walkers
   1    continue
        nac=nac+1
        iclx(nac)=iatx
        icly(nac)=iaty
        if(nac<nacf) then ! execute until reach final number in cluster
          goto 2
        endif
        ! dimensionality check
        irm=50
        do i=2,irm ! radius from center, square
          is=0 ! sum of mass
          rd=i ! radius
          do j=1,nacf
            x=iclx(j) ! loop over atoms in cluster
            y=icly(j)
            if((abs(x)<rd).and.(abs(y)<rd)) then
              is=is+1 ! add atom to sum
            endif
          enddo
          sm=is
          write(8,112) log(rd),log(sm)
        enddo
        do j=1,nacf
          write(7,111) iclx(j),icly(j)
        enddo
        end



