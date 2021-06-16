        program laplace ! will be simulating capacitor
        real v(100,100),v1(100,100) ! v is initial v1 next
        print*, "Enter the number of iterations"
        read(5,*) ni !# of iterations 
        v=0.0 !array value cast
        v1=0.0
        open(7,file="potential")
  111   format(100(2x,e12.5))!40 lines, 40 elements each
        do j=40,60 ! initial conditions
          v(40,j)=-1.0
          v1(40,j)=-1.0
          v(60,j)=1.0
          v1(60,j)=1.0
        enddo
        ! store
        do ii=1,ni ! iterations, main run
          call update(v,v1,100) ! v1 is output
          do j=40,60
            v1(40,j)=-1.0
            v1(60,j)=1.0
          enddo
          call update(v1,v,100) ! v is output
          do j=40,60
            v(40,j)=-1.0
            v(60,j)=1.0
          enddo
          dv = 0.0 !maximum difference between prev and next iteration
          do i=2,99 !cover all i,j
            do j=2,99
              err=abs(v1(i,j)-v(i,j))
              if(err>dv)then
                dv=err
              endif
            enddo
          enddo
          if(dv<10.0e-5) then !converged
            goto 1
          endif
        enddo ! iterations
  1     continue ! target for exit condition
        if(dv>10.0e-5) then
          print*, "Potential has not been converged"
          print*, "dv=",dv
        endif
        do i=1,100
          write(7,111) (v(i,j),j=1,100)
        enddo
        end

        subroutine update(v,v1,n) ! input pot, output pot, #points
        real v(n,n),v1(n,n)
        do i=2,n-1 ! 2 because iteration steps
          do j=2,n-1
            v1(i,j)=0.25*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1))
          enddo
        enddo
        return
        end



