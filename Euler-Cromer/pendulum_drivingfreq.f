        program pendulum
        real ant(300), omega(300), theta(300), t(300) ! ant is analytical theta
        open(7,file="pen_out_Cromer_friction_drivingfreq")
  111   format(3(2x,e12.5))
        print*, "Enter the friction coefficient"
        read(5,*) q
        print*, "Enter the driving frequency value"
        read(5,*) df
        ant=0.0 ! all elements of array set to 0, just in case
        omega=0.0
        theta0=3.0 ! 3 degrees, initial angle
        pi=3.1415927
        theta0=theta0*pi/180.0 ! convert to radians
        l=1.0 ! length
        g=9.8
        freq=sqrt(g) ! natural frequency, l=1, not included
        period=2.0*pi/freq
        ti=6.0*period ! time interval, 6 periods considered
        t(1)=0.0
        theta(1)=theta0
        dt=ti/300.0
        omega(1)=0.0
        t(1)=0
        do i=1,299 ! it=timestep, last will end up being 300
          t(i+1)=t(i)+dt
          ant(i+1)=theta0*cos(freq*t(i+1))
          omega(i+1)=omega(i)-(g*theta(i)+q*omega(i)-
     *    df*cos(freq*t(i+1)))*dt ! freq is natural freq
          theta(i+1)=theta(i)+omega(i+1)*dt !omega(i)->omega(i+1)Cromer
          write(7,111) t(i+1),ant(i+1),theta(i+1)
        enddo
        end
