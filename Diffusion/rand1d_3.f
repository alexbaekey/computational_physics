        program rand1d ! add many workers"
        print*, "Enter the number of steps"
        read(5,*) ns
        print*, "Enter the number of walkers"
        read(5,*) nw
        open(7,file="walk_2")
  111   format(2(2x,i4)) ! 4 available digits
        do iw=1,nw ! walkers
          iy=0.0 ! initial position
          do i=1,ns
            r1=rand(0) ! random number generator, 0 means new every time
            if(r1<0.5) then
              iy=iy+1
            else
              iy=iy-1
            endif
            write(7,111) i,iy ! step number, displacement
        enddo
        enddo ! walkers
        end
