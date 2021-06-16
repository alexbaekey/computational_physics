        program eden ! no diffusion
        integer iclx(1000),icly(1000),ipx(1000),ipy(1000)
        open(7,file="cluster")
        open(8,file="perimeter")
        open(9,file="dimf") 
  111   format(2(2x,i4)) ! 2 clolumn, x and y, 2 spaces, 4 digits
  112   format(2(2x,e12.5)) 
        print*,"Enter the cluster size"
        read(5,*) nacf
        iclx(1)=100
        icly(1)=100
        nac=1 ! number of cluster atoms
   1    continue
        nap=1 ! number of perimeter atoms
        do i=1,200 !test positions around x
          do j=1,200 !around y
            is=0
            do ic=1,nac ! check whether already occupied
              if((i==iclx(ic)).and.(j==icly(ic))) then
                is=is+1
              endif
            enddo
            if(is==0) then !decision about perimeter position
              do ic=1,nac ! check if neighboring
                ix=abs(iclx(ic)-i)
                iy=abs(icly(ic)-j)
                if(((ix==1).and.(iy==0)).or.((ix==0).and.(iy==1))) then
                  ipx(nap)=i
                  ipy(nap)=j
                  nap=nap+1
                endif
              enddo
            endif
          enddo ! end around y
        enddo ! end around x
        r1=rand(0) ! generate random number
        ic=int(r1*(nap-1)) ! site to be occupied
        if(ic==0) then
          ic=1
        endif
        nac=nac+1 ! adding atom to the cluster
        iclx(nac)=ipx(ic)
        icly(nac)=ipy(ic)
        if(nac<nacf) then
          goto 1
        endif
        do i=1,nacf
          write(7,111) iclx(i),icly(i)
        enddo
        do i=1,nap-1 !ignore extra nap that wasn't assigned at end
          write(8,111) ipx(i),ipy(i)
        enddo
        ! dimensionality check
        irm=50
        do i=2,irm ! radius from center, square
          is=0 ! sum of mass
          rd=i ! radius
          do j=1,nacf
            x=iclx(j) ! loop over atoms in cluster
            y=icly(j) 
            if((abs(x)<rd+100).and.(abs(y)<rd+100)) then
              is=is+1 ! add atom to sum
            endif
          enddo
          sm=is
          write(9,112) log(rd),log(sm)
        enddo 
        end




