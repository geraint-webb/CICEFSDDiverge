!GW
module NilsFSD

  use icepack_kinds
      use icepack_parameters, only: p01, p5, c0, c1, c2, c3, c4, c10
      use icepack_parameters, only: bignum, puny, gravit, pi
      use icepack_warnings, only: warnstr, icepack_warnings_add,  icepack_warnings_aborted


      implicit none
      public :: correct_FSD


    contains

      subroutine correct_FSD(ncat, nfsd, trcrn, aicen, divu, floe_rad_c, floe_binwidth, floe_rad_l, dt, hin_max,d_afsd_nils)

   integer (kind=int_kind), intent(in) :: ncat, nfsd
    real (kind=dbl_kind), intent(in) :: divu, dt
    integer (kind=int_kind) :: counter, n, nnum, k, closest_cat, p, points, max_checks, l, integrate_points, d, last_cat_hold
    real (kind=dbl_kind) :: ice_remaining, tol, factor, divu_per_day, tester, icelost, weight, time_steps_per_day, max_cat_value, floe_diameter, floe_increment, a0_search, d0_search_last, dP_search, dP_search_deriv, min_bnd, max_bnd, diam, dc, di, point1, point2, ai, ac, dP_last, grid_size
    real (kind=dbl_kind), dimension(:,:), intent(inout) :: trcrn           ! tracer array - just fsd
    real (kind=dbl_kind), dimension(:), intent(in) :: aicen, floe_rad_c, floe_binwidth, floe_rad_l, hin_max
    real (kind=dbl_kind), dimension(ncat-1) :: dhmin_max
    integer (kind=int_kind), dimension(ncat) :: last_cat
    real (kind=dbl_kind), dimension(3) :: aCof, eCof, bCof, dCof
    real (kind=dbl_kind), dimension(nfsd) :: dP, dPnorm, new_dP, dP_afsd
    real (kind=dbl_kind), dimension(ncat) :: run_ice_cat, factor_change
    real (kind=dbl_kind), dimension(nfsd) :: closest_boundary, afsd_weight, diff
    real (kind=dbl_kind), dimension(nfsd+1) :: fsd_bounds, a_bounds
    real (kind=dbl_kind), dimension(nfsd,ncat) :: final_afsd, afsd, d_afsdn_nils, dP_saved, differ, last_cat_location, dP_afsd_kn, int_afsd
    real (kind=dbl_kind), dimension(nfsd), intent(inout) :: d_afsd_nils
    real (kind=dbl_kind) :: a, e, b, d0, dPmean, run_ice_total
    logical, dimension(ncat) :: cat_flag_found, saturated
    logical :: flag_warning, max_checks_reached
    logical, dimension(nfsd) :: empty_cats
      max_checks_reached=.false.
      flag_warning = .false.
      d_afsdn_nils=c0
      int_afsd=trcrn(:,:)
      d_afsd_nils=c0
      divu_per_day=divu*(86400)
      if (divu_per_day>0.055) then
         divu_per_day=0.0549
      else if (divu_per_day>0.003 .and. divu_per_day<0.009) then
         divu_per_day=0.003
      end if
    
      if (divu_per_day>0 .and. divu_per_day<=0.055) then ! Theory only holds for this region... 
      bCof=(/-1.27e-1,7.56e-3,-2.74e-5/)
      aCof=(/4.72e6,3.3e3,-3.89e2/)
      eCof=(/6.59e2,-3.87e1,-5e-1/)
      dCof=(/1.77e6,-1.19e5,3.98e4/)

      a=aCof(1)*divu_per_day**2 +aCof(2)*divu_per_day +aCof(3)
      e=eCof(1)*divu_per_day**2 +eCof(2)*divu_per_day +eCof(3)
      b=bCof(1)*divu_per_day**2 +bCof(2)*divu_per_day +bCof(3)

      time_steps_per_day=86400/dt! Only want to apply once per day
      ! Calculate d0
      d0_search_last=c0

      dP_search=c0
      max_checks=2000

      !Forcing bounds
      fsd_bounds(1:nfsd)=2*floe_rad_l(:)
      if (2*floe_rad_l(nfsd)<50000) then
         fsd_bounds(nfsd+1)=50000
      else
         fsd_bounds(nfsd+1)=2*floe_rad_l(nfsd)! Arbitrary rework for 24 bins
      end if        
    
      ! Convert to area
      a_bounds(:)=0.66*fsd_bounds**2

  
      min_bnd=a_bounds(1)
      max_bnd=a_bounds(nfsd+1)

      a0_search=(max_bnd+min_bnd)/2

      ! Simple bisection method
      dP_search=c0
      tol=10**4 ! Tolerance for convergence relativ eto the amount on order e8
      if (fsd_bounds(nfsd+1) .NE. 50000) then
      integrate_points=200
      diam=(max_bnd-min_bnd)/integrate_points
      ai=a_bounds(1)
      do p=1, max_checks
         dP_search=c0
         ac=c0
         do d=2, integrate_points
               if (d==2) then
                  point1=a*(sqrt((1/0.66)*ai)**e)*tanh(b*(sqrt((1/0.66)*a0_search)-sqrt((1/0.66)*ai)))*0.66*(sqrt((1/0.66)*ai))**2
               end if
               ac=ac+diam
               point2=a*(sqrt((1/0.66)*ac)**e)*tanh(b*(sqrt((1/0.66)*a0_search)-sqrt((1/0.66)*ac)))*0.66*sqrt((1/0.66)*ac)**2
               dP_search=dP_search+(point1+point2)*diam/2
               point1=point2
         end do
         if (abs(dP_search) <= tol) then
            exit
         else if (dP_search<0) then
            min_bnd=a0_search
            a0_search=(max_bnd+min_bnd)/2
         else if (dP_search > 0) then
            max_bnd=a0_search
            a0_search=(max_bnd+min_bnd)/2
         else if (p==max_checks .and. abs(dP_search) > tol) then 
            write(warnstr,*) dP_search, 'Nils d0 failed to converge, GW', a0_search, divu_per_day
            call icepack_warnings_add(warnstr)
         end if
      end do
   else
      ! For 12 and 16 cats just use Hutter's fit
      a0_search=0.66*(dCof(1)*divu_per_day**2 +dCof(2)*divu_per_day +dCof(3))**2
      end if

      ! Error message if d0 doesn't converge
      if (dP_search >tol .or. dP_search<-tol) then
         write(warnstr,*) dP_search, 'Nils d0 failed to converge, GW', a0_search, divu_per_day
         call icepack_warnings_add(warnstr)
         a0_search=975062402.2737021! Temp fix preventing crash
      end if
      




      ! Now we have d0 we can calculate the mean dP for each bin using the trapezoidal rule
      dP=c0
      grid_size=200
      do k=1, nfsd-1 ! For each bin
      diam=real((a_bounds(k+1)- a_bounds(k))/grid_size,8)
      ai=a_bounds(k) ! Start at the lower bound
      ai=real(ai,8)
      ac=ai
      
      dP_search=c0
      do d=2, grid_size ! For each point in the bin
         if (d==2) then ! Only calculate the first point once
            point1=(a*(sqrt((1/0.66)*ai)**e)*tanh(b*(sqrt((1/0.66)*a0_search)-sqrt((1/0.66)*ai)))*0.66*(sqrt((1/0.66)*ai))**2) 
            point1=real(point1, kind=8)/(real(5e7, kind=8))! Converts from per km^2 to per m^2 
            point1=point1/time_steps_per_day ! Converts from per day to per timestep
               end if
               ac=real(ac+diam, kind=8) ! Move to the next point
               point2=(a*(sqrt((1/0.66)*ac)**e)*tanh(b*(sqrt((1/0.66)*a0_search)-sqrt((1/0.66)*ac)))*0.66*sqrt((1/0.66)*ac)**2)
               point2=real(point2, kind=8)/(real(5e7, kind=8))
               point2=point2/time_steps_per_day
               dP_search=dP_search+(point1+point2)*diam/2 ! Implement trapezoidal rule
               ! Error messages catching potential mistakes
               if (point1<0 .or. point2<0) then
                  write(warnstr,*) 'point1/point2 less than 0, GW', k, d, point1, point2, tanh(b*(sqrt((1/0.66)*a0_search)-sqrt((1/0.66)*ac))), divu_per_day ! The other way to fix this error is instead of setting intercept to zero  once we get Nil's coefficents just find the min divergencee for this rule to hold
                  call icepack_warnings_add(warnstr)
               end if
               if (ac>a_bounds(k+1)) then
                  write(warnstr,*) 'ac>a_bounds(k+1)', k, d                                                                                                                                                     
                  call icepack_warnings_add(warnstr)
               end if 
               point1=point2
               
            end do
            dP(k)=dP_search/(a_bounds(k+1)- a_bounds(k)) ! Dividing by delta a to get the mean value
      
            if (dP(k)<0 .or. abs(dP(k))>1) then ! Sanity Check
               write(warnstr,*) 'dP(k)<0 or dP(k)>1', dP(k), dP_search, k, divu_per_day, a0_search, time_steps_per_day
               call icepack_warnings_add(warnstr)
            end if
         end do  
       
         dP_last=-sum(real(dP,8)) ! Sum up the mean values to get the total negative change

         if (dP_last>0 .or. abs(dP_last)>1) then ! Sanity Check
            write(warnstr,*) 'dP_last', dP_last, divu_per_day, dP(nfsd-1),dP(nfsd-2),dP(1), a0_search, diam, time_steps_per_day, point1, a
            call icepack_warnings_add(warnstr)
         end if
         ! Update dP to include the last bin
       dP(nfsd)=dP_last
      
       ! Find last category for each ncat
       cat_flag_found=.true. ! Setting all of array to true
       run_ice_total=c0
       dP_saved=c0
       do n=1, ncat ! Setting each categories FSD to the initial FSD - note these categories will change at different rates
          dP_saved(:,n)=dP
       end do   
       counter=c0
       run_ice_cat=c0
       differ=c0
       ! #bins (nfsd) x #cats (ncat)
       
      saturated=.FALSE. ! Keeping track of which categories there is potential ice left to move
      dP_afsd_kn=c0
       do while (run_ice_total<=-0.99*ncat*dP_last .and. any(cat_flag_found) .and. counter<=nfsd .and. any(.not. saturated))
       last_cat_hold=c0
       last_cat=c0
       cat_flag_found=.false.
       do n=1, ncat
          last_cat_hold=1
          do k=1, nfsd
             if (a_bounds(k)> 4*0.66*100**2 .and. trcrn(k,n)>10**-2) then 
                last_cat_hold=k
                cat_flag_found(n)=.true.
                end if
             end do
          last_cat(n)=last_cat_hold
       end do
       
       do n=1, ncat
          dP_afsd(:)=c0
          
          if (cat_flag_found(n)==.true. .and. saturated(n)==.false.) then ! If not true once add left over ice to be moved for this category into run_ice_total to end the loop... this should fix the bhug 
             new_dP=c0
          do k=1, nfsd! Save intial histogram but force it to not use the last?
             if (k<last_cat(n)) then   
                new_dP(k)=dP_saved(k,n)
             else if (k==last_cat(n)) then
                new_dP(k)=sum(dP_saved(last_cat(n):,n))
                if (new_dP(k)>0) then
                   write(warnstr,*) 'last dP(k)>0!',new_dP(k)
         call icepack_warnings_add(warnstr)
                end if
             end if
          end do
          weight=c0
          icelost=c0
          do k=1, nfsd
             if (new_dP(k)<0 .and. abs(new_dP(k))>trcrn(k,n)) then
                 icelost=icelost+abs(trcrn(k,n))                                                        
                 dP_afsd(k)=-abs(trcrn(k,n))
                 d_afsdn_nils(k,n)=d_afsdn_nils(k,n) -abs(trcrn(k,n))
             else if (new_dP(k)<0 .and. abs(new_dP(k))<=trcrn(k,n)) then
                icelost=icelost+abs(new_dP(k))                   ! Should exit loop
                saturated(n)=.true.
                dP_afsd(k)=new_dP(k)
                d_afsdn_nils(k,n)=d_afsdn_nils(k,n)+new_dP(k)
             else
                weight=weight+new_dP(k)
                if (new_dP(k)<0) then
                   write(warnstr,*) 'new_dP(k)<0 when should be greater, GW', new_dP(k), dP_search
                   call icepack_warnings_add(warnstr)
                end if
             end if
          end do
          ! Now distrubute
          do k=1, nfsd
             if (new_dP(k)>0) then
                dP_afsd(k)=(new_dP(k)/weight)*icelost
                d_afsdn_nils(k,n)=d_afsdn_nils(k,n)+(new_dP(k)/weight)*icelost
             end if
          end do
          dP_afsd_kn(:,n)=dP_afsd_kn(:,n)+dP_afsd(:)
          if (last_cat(n) .ne. nfsd) then
            run_ice_cat(n)=run_ice_cat(n)+icelost+(icelost/abs(new_dP(last_cat(n))))*sum(dP_saved(last_cat(n):ubound(dP_saved, 1)-1,n))!(icelost/-new_dP(last_cat(n)))
          else
            run_ice_cat(n)=run_ice_cat(n)+icelost
          end if

          if (icelost<0) then
               write(warnstr,*) 'icelost<0, GW', icelost, n, counter
               call icepack_warnings_add(warnstr)
          end if

          if (sum(trcrn(:,n))>1.03*sum(trcrn(:,n)+dP_afsd(:)) .or. sum(trcrn(:,n))<0.97*sum(trcrn(:,n)+dP_afsd(:))) then
             write(warnstr,*) 'area not conserved, GW', sum(trcrn(:,n)), sum(trcrn(:,n)+dP_afsd(:)), counter
                   call icepack_warnings_add(warnstr)
          end if
        
         trcrn(:,n)=trcrn(:,n)+dP_afsd(:)  
         factor_change(n)=1-(run_ice_cat(n)/(abs(dP_last)))
         if (factor_change(n)<0 .and. saturated(n)==.false.) then
         write(warnstr,*) 'Factor Change negative, GW', factor_change, n, counter
         call icepack_warnings_add(warnstr)
         end if
         if (run_ice_cat(n)>1.01*abs(dP_last)) then
            write(warnstr,*) 'Too mucg ice, GW', run_ice_cat(n), -dP_last, n, counter
            call icepack_warnings_add(warnstr)        
         end if
         dP_saved(:,n)=dP*factor_change(n)
     
    else
       run_ice_cat(n)=abs(dP_last)
       saturated(n)=.true.
    end if             
    end do
    
    run_ice_total=sum(run_ice_cat)
   !if (counter>0) then
   !      write(warnstr,*) 'GW', counter, sum(run_ice_cat), -0.99*ncat*dP_last, factor_change(ncat), factor_change(1)
   !      call icepack_warnings_add(warnstr)
   !   end if

    counter=counter+1
    if (counter==nfsd+1) then
       write(warnstr,*) 'Counter =13!, GW', dP_saved(nfsd,ncat),run_ice_total, -0.9*ncat*dP_last, factor_change(ncat), factor_change(1)
       call icepack_warnings_add(warnstr)
    end if
   
 end do

 ! This is one of the last test cases I have setup to check if the code is working correctly
 if (run_ice_total<-0.9*ncat*dP_last .and. any(.not. saturated)) then
    write(warnstr,*) 'run_ice_total<-0.9*ncat*dP_last, GW', run_ice_total, -0.9*ncat*dP_last
    call icepack_warnings_add(warnstr)
 end if
 

 do n=1, ncat
    if (sum(trcrn(:,n))>puny .and. (sum(trcrn(:,n))>1.03*sum(int_afsd(:,n)) .or. sum(trcrn(:,n))<0.97*sum(int_afsd(:,n)))) then
             write(warnstr,*) 'area not conserved final check, GW', n, sum(trcrn(:,n)), sum(int_afsd(:,n))
                   call icepack_warnings_add(warnstr)
          end if
 end do

 do n=1, ncat
 if (sum(trcrn(:,n))>1.01) then
          write(warnstr,*) 'Ice area not conserved, GW', n, sum(trcrn(:,n)), sum(int_afsd(:,n))
          call icepack_warnings_add(warnstr)
       end if
end do 
 ! Archive the final FSD
d_afsdn_nils=c0
d_afsdn_nils=trcrn-int_afsd
do k = 1, nfsd
         d_afsd_nils(k) = c0
         do n = 1, ncat
            d_afsd_nils(k) = d_afsd_nils(k) + aicen(n)*d_afsdn_nils(k,n)
         end do ! n                                                                                                                                                                                                                              
      end do 
    end if
 
end subroutine

end module


    
