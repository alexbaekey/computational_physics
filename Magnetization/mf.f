        program meanfield
        real smag(200)
        dt=0.04 ! temperature step
        eps=1.0e-4
        open(7,file="mag-mf")
  111   format (2(2x,e12.5))
        do it=1,200 ! loop over temps
          sm=0.5 ! guess the spin, sm is ave mag moment
          rt=it*dt ! temperature in terms of J/k_B
   1      continue
          rtnh=tanh(4.0*sm/rt)
          fs=sm-rtnh
          dfs=1.0-4.0*(1.0-rtnh**2)/rt
          if(dfs==0) then
            dfs=1.0e-7
          endif
          sm0=sm !sm0 is s_i, sm is s_i+1
          sm=sm0-fs/dfs
          delta=abs(sm-sm0)
          if(delta>eps) then
            goto 1
          endif
          smag(it)=sm
          write(7,111) rt,smag(it)
        enddo
        end
