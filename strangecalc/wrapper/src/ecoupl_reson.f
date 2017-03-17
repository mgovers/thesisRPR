

c	FORTRAN CODE TO COMPUTE EM FORMFACTORS OF RESONANCES
c	
c	OBTAINED FROM https://userweb.jlab.org/~isupov/couplings/



c###############################################################################################


        

c       example of using ecouple function
c       calculating A32 of D33(1700) at Q2 = 2.1 GeV2 

c       program reson
c	implicit none
c	REAL Q2,result,ecoupl
c	result =  ecoupl(2,13,2.1)   
c	print *,result
c	end



c###############################################################################################


c        The function calculates electrocouplings of N*
c		 as a function of Q2 using phenomenological fits.

c        The first argument is the helecity amplitude
c        i=1  corresponds to A12
c        i=2  corresponds to A32
c        i=3  corresponds to S12  

c        The second argument is the resonance number
c        25  is   p33(1232)
c        20  is   p11(1440)
c	 15  is   d13(1520)
c        18  is   s11(1535)
c	 21  is   s31(1620)  	 
c	 19  is   s11(1650) 
c        6   is   f15(1685)
c        13  is   d33(1700)
c        17  is   p13(1720)
c        14  is   p13(1720) missing state

c		The third argument is Q2 from 0 to 5.0 GeV2  
  
        
        real function ecoupl(i,j,Q2)  ! i=1 -A12,  i=2 -A32, i=3 -S12  j-resonance number
        implicit none
        real Q2,result,value,A12,A32,S12
        integer i,j
        
        ecoupl = 0
        value = 0
        
        if(Q2.lt.0.or.Q2.gt.5.0) return
        	
        
        if(j.eq.14) then  ! missing p13(1720)
        if(i.eq.1) value =  (36.5+1044.05*Q2)/(1.+19.5*Q2*sqrt(Q2)) 
c		if(i.eq.2) value = (-39.6+101.4*Q2-119.1*Q2*Q2)/(1.+1.37*Q2**4)
		if(i.eq.2) value = (-37.9)/(1.+0.455*Q2**2*sqrt(Q2))
		if(i.eq.3.and.Q2.gt.0.) value =  (93.056)/(1.+1.61*Q2*sqrt(Q2)) 
        ecoupl = value     
        return      
        endif
        
        if(j.eq.17) then  !  p13(1720) 
        if(i.eq.1) value =  (90)/(1.+3.16*Q2*sqrt(Q2)) 
		if(i.eq.2) value = (-35.87-6.85*Q2)/(1.+ 0.118*Q2**4)
        if(i.eq.3.and.Q2.gt.0.) value=-3.09/(1.-3.7*Q2+2.86*Q2*sqrt(Q2))
        ecoupl = value     
        return      
        endif
        
        if(j.eq.13) then  !  d33(1700) 
        if(i.eq.1) value =  (118.28)/(1.+0.72*Q2*Q2+3.26*Q2**2.5) 
		if(i.eq.2) value = (101.)/(1.+16.6*Q2*Q2-5.0*Q2**2.5)
		if(i.eq.3.and.Q2.gt.0.) value =  (24.5)/(1.+0.15*Q2*sqrt(Q2)) 
        ecoupl = value     
        return      
        endif
        
        if(j.eq.6) then  !  f15(1680) 
        
        A12=1.+3.04*Q2+0.0305*Q2**2
        A12=A12/(1.+0.034*(Q2-2.)**2)
        A12=A12/(1.+0.765*Q2*sqrt(Q2))
        A12=-20.*A12/(1.+0.09*Q2**2)
        if(i.eq.1) value = A12
        
        A32=4.83*Q2-2.68*Q2**2+1.1*Q2**2*sqrt(Q2)
        A32=134.2/(1.+A32)
        if(i.eq.2) value = A32
        
        S12=-109.04*(1.+Q2)
        S12=S12/(1.+1.9*Q2**4)
		if(i.eq.3.and.Q2.gt.0.) value =  S12
		
        ecoupl = value     	    
        return
        endif
        
        if(j.eq.19) then  !  s11(1650) 
c        if(i.eq.1) value =  (551.07)/(1.+17.73*Q2*Q2) 
		if(i.eq.1) value=(47.4-19.6*Q2)/(1.-1.46*Q2*sqrt(Q2)+1.17*Q2**3)
		if(i.eq.3.and.Q2.gt.0.) value=-2.67/(1.-2.82*Q2+2.0*Q2*sqrt(Q2)) 
        ecoupl = value     
        return      
        endif
        
        if(j.eq.21) then  !  s31(1620) 
        if(i.eq.1) value =  (47.2)/(1.+3.71*Q2*sqrt(Q2))
		if(i.eq.3.and.Q2.gt.0.) value =  (-25.)/(Q2*sqrt(Q2)) 
        ecoupl = value     
        return      
        endif
        
     	if(j.eq.18) then  !  s11(1535) 
     	
     	A12=92.5029+1.45023*Q2
        A12=A12/(1.+0.1095*Q2**2-0.000322*Q2**2*sqrt(Q2))
        if(i.eq.1) value = A12        
        
        S12=-9.758811-4.231412*Q2
        S12=S12/(1.-0.7341952*Q2**2+0.5087887*Q2**2*sqrt(Q2))
		if(i.eq.3.and.Q2.gt.0.) value =  S12
		        
        ecoupl = value              
        return
        endif
        
        
     	if(j.eq.15) then  !  d13(1520) 
        
        A12= -23.357-151.199533*Q2
        A12=A12/(1.+2.01489898*Q2**2-0.2654327*Q2**2*sqrt(Q2))
        A12=A12*0.9
        if(i.eq.1) value = A12
        
        A32=3.322979*Q2-2.0339966*Q2**2+1.622563*Q2**2*sqrt(Q2)
        A32=162.458285/(1.+A32)
        A32=A32*0.9
        if(i.eq.2) value = A32
        
        S12=1.73*Q2-2.8*Q2**2+2.91*Q2**2*sqrt(Q2)
        S12=-67.32/(1.+S12)
		if(i.eq.3.and.Q2.gt.0.) value =  S12        
                                        
        ecoupl = value              
        return
        endif
        
     	if(j.eq.20) then  !  p11(1440)
        
        A12=-68.7866+21.3966*Q2+79.8415*sqrt(Q2)
        A12=A12/(1.-0.7178*Q2**2+0.5663*Q2**2*sqrt(Q2))
        if(i.eq.1) value = A12
        
        S12=31.19227+3.53338*Q2
        S12=S12/(1.-0.278265*Q2**2+0.3677575*Q2**2*sqrt(Q2))
		if(i.eq.3.and.Q2.gt.0.) value =  S12                      	
     	         
        ecoupl = value              
        return
        endif
        


	if(j.eq.25) then  !  p33(1232) 

        A12= -170.06/((1.+Q2)*(1.+0.1609*Q2**2-0.002*Q2**4))
        if(i.eq.1) value = A12
        

        
        A32= -321.06/((1.+Q2)*(1.+0.16*Q2**2-0.002*Q2**4))
        if(i.eq.2) value = A32
        

        
        S12= 29.76/((1.+Q2)*(1.+0.0135*Q2**2-0.00046*Q2**4))
		if(i.eq.3.and.Q2.gt.0.) value =  S12        
                                        
        ecoupl = value              
        return
        endif        
        
        
        
        
             
             	


        ecoupl = value              
        return        
        end
