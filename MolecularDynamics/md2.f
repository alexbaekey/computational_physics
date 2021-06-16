        program md
        real x(16),y(16),vx(16),vy(16),frcx(16),frcy(16),xch(3,16),
     *  ych(3,16),dr(16,16),rdf(150) ! xch and ych are xchange, ychange,
c       dr is distance between all atoms, rdf is radial distrib function
        open(7,file="xy",status="replace")
  111   format(2(2x,e12.5))
        open(8,file="coords",position="rewind")
        open(9,file="raddf",status="replace")
        print*, "Enter v0"
        read(5,*) v0 ! maximum speed
        per=4.0 ! periodicity
        dt=0.002 ! timestep
        print*, "Enter the number of time steps"
        read(5,*) itn
        print*,"Enter the time interval for saving"
        read(5,*) ibeg,iend
        print*, "How to start: scratch=0/coords=1"
        read(5,*) istart
! initialization
        if(istart==0) then
          print*, "Start from srcatch"
          k=0
          do i=1,4
            do j=1,4
              k=k+1
              x(k)=0.5+j-1
              y(k)=0.5+i-1
            enddo
          enddo
        else
          print*, "Start from file"
          do i=1,16
            read(8,111) x(i),y(i)
          enddo
        endif ! start from scratch
        do ia=1,16
          r3=0.5-rand(0)
          vx(ia)=1.414214*v0*r3
          r4=0.5-rand(0)
          vy(ia)=1.414214*v0*r4
        enddo
        svx=0.0
        svy=0.0
        do ia=1,16
          svx=svx+vx(ia)
          svy=svy+vy(ia)
        enddo
        print*,svx,svy
        vx16=svx/16.0
        vy16=svy/16.0
        do ia=1,16
          vx(ia)=vx(ia)-vx16
          vy(ia)=vy(ia)-vy16
        enddo
        svx=0.0
        svy=0.0
        do ia=1,16
          svx=svx+vx(ia)
          svy=svy+vy(ia)
        enddo
        do ia=1,16
          xch(2,ia)=x(ia)
          ych(2,ia)=y(ia)
          xch(1,ia)=xch(2,ia)-vx(ia)*dt
          ych(1,ia)=ych(2,ia)-vy(ia)*dt
        enddo
c end of initialization
        rdf=0.0
        do it=1,itn !!! time steps
          call forces(x,y,frcx,frcy,16,per,dr)
          do ia=1,16 ! atoms   
            xch(3,ia)=2.0*xch(2,ia)-xch(1,ia)+frcx(ia)*dt*dt
            ych(3,ia)=2.0*ych(2,ia)-ych(1,ia)+frcy(ia)*dt*dt
            x(ia)=xch(3,ia)
            y(ia)=ych(3,ia)
c now shift
            xch(1,ia)=xch(2,ia)
            ych(1,ia)=ych(2,ia)
            xch(2,ia)=xch(3,ia)
            ych(2,ia)=ych(3,ia)
c periodic boundary condition 
            if(x(ia)>per) then
              x(ia)=x(ia)-per
              xch(1,ia)=xch(1,ia)-per
              xch(2,ia)=xch(2,ia)-per
            else
              if(x(ia)<0) then
                x(ia)=x(ia)+per
                xch(1,ia)=xch(1,ia)+per
                xch(2,ia)=xch(2,ia)+per
              endif
            endif
            if(y(ia)>per) then
              y(ia)=y(ia)-per
              ych(1,ia)=ych(1,ia)-per
              ych(2,ia)=ych(2,ia)-per
            else
              if(y(ia)<0) then
                y(ia)=y(ia)+per
                ych(1,ia)=ych(1,ia)+per
                ych(2,ia)=ych(2,ia)+per
              endif
            endif
          enddo !! atoms
          if((it>ibeg).and.(it<iend)) then
            do irdf=1,150 ! loop over radius
              do i=1,16
                ra=irdf*0.02
                dr7=dr(7,i)
                if((dr7>ra).and.(dr7<=ra+0.02)) then
                  rdf(irdf)=rdf(irdf)+1.0
                endif
              enddo
            enddo ! end radius loop
            if(mod(it,10)<1) then
              do i=1,16
                write(7,111) x(i),y(i)
              enddo
            endif
          endif
        enddo !!time steps
        do i=1,16
          write(8,111) x(i),y(i)
        enddo
        do i=1,150
          write(9,111) i*0.02,rdf(i)
        enddo
        end

        subroutine forces(x,y,frcx,frcy,n,per,dr)
        real x(n),y(n),frcx(n),frcy(n),dr(n,n),dx(n,n),dy(n,n)
        do i=1,n
          do j=i,n ! translations
            dxc=x(i)-x(j)
            dyc=y(i)-y(j)
            ddx=dxc
            ddy=dyc
            drc=sqrt(dxc**2+dyc**2) !without translation
            dd=sqrt((ddx-per)**2+ddy**2) ! dx-L
            if(dd<drc) then
              drc=dd ! keep smallest seperation
              dxc=ddx-per
              dd1=sqrt(dxc**2+(ddy-per)**2)
              if(dd1<drc) then
                drc=dd1
                dyc=ddy-per
              endif
              dd2=sqrt(dxc**2+(ddy+per)**2)
              if(dd2<drc) then
                drc=dd2
                dyc=ddy+per
              endif
            endif ! dx-L
            dd=sqrt((ddx+per)**2+ddy**2) ! dx+L
            if(dd<drc) then
              drc=dd
              dxc=ddx+per
              dd1=sqrt(dxc**2+(ddy-per)**2) 
              if(dd1<drc) then
                drc=dd1
                dyc=ddy-per
              endif
              dd2=sqrt(dxc**2+(ddy+per)**2)
              if(dd2<drc) then
                drc=dd2
                dyc=ddy+per
              endif
            endif !!dx+L
            dd=sqrt(ddx**2+(ddy-per)**2)
            if(dd<drc) then
              drc=dd
              dxc=ddx
              dyc=ddy-per
            endif            
            dd=sqrt(ddx**2+(ddy+per)**2)
            if(dd<drc) then
              drc=dd
              dxc=ddx
              dyc=ddy+per
            endif
            dr(i,j)=drc
            dx(i,j)=dxc
            dy(i,j)=dyc
          enddo
        enddo
        do i=1,n
          do j=1,i-1
            dr(i,j)=dr(j,i)
            dx(i,j)=-dx(j,i)
            dy(i,j)=-dy(j,i)
          enddo
        enddo
        do i=1,n
          sfx=0.0
          sfy=0.0
          do j=1,n
            if(i.ne.j) then ! not equal operator
              fm=24.0*(2.0/dr(i,j)**13-1.0/dr(i,j)**7)
              sfx=sfx+fm*dx(i,j)/dr(i,j)
              sfy=sfy+fm*dy(i,j)/dr(i,j)
            endif
          enddo
          frcx(i)=sfx
          frcy(i)=sfy
        enddo    
        return
        end





