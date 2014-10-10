SUBROUTINE Decisions

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

! Local Variables
integer ::  it,ia,ip,ie,iz,ie2,iz2,BLbind,imnext(2),im
REAL(8)			::  emuc(ngpa),lpmnext(2),lnextm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
! Start with last period, eat everything  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			

if (Display==1) write(*,*) 'Solving for decision rules at age ', Ttot
			
DO ip = 1,ngpp
    conret(Tret,:,ip) = xgridret(Tret,:,ip)
    IF (QuadraticPref==1) THEN
    	mucret(Tret,:,ip) = qprefa*(qprefb - conret(Tret,:,ip))			
    ELSE
    	mucret(Tret,:,ip) = conret(Tret,:,ip)**(-gam)			
    END IF
END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retired Agents: Pension value as state variable !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO it = Tret-1,1,-1
	if (Display==1) write(*,*) 'Solving for decision rules at age ', Twork+it				

    !$OMP PARALLEL DO PRIVATE(BLbind)
	DO ip = 1,ngpp				
		!solve on tt+1 grid
		mucret1(it,:,ip) = bet*(surprob(it)/annprem(it))*Rnet*mucret(it+1,:,ip)		
		IF (QuadraticPref==1) THEN
		    conret1(it,:,ip) = qprefb - mucret1(it,:,ip)/qprefa
		ELSE
		    conret1(it,:,ip) = mucret1(it,:,ip)**(-1.0/gam)  
		END IF
		assret1(it,:,ip) = (conret1(it,:,ip) +agridret(it+1,:,ip)*annprem(it)-pgrid(it,ip) -tranret(it,:,ip))/Rnet

		!deal with borrowing limits
		IF (MINVAL(agridret(it,:,ip)) >= assret1(it,1,ip)) THEN		!BL does not bind anywhere
			BLbind = 0
		ELSE				!find point in tt grid where BL starts to bind
			BLbind = MAXLOC(agridret(it,:,ip), DIM=1,MASK = agridret(it,:,ip) .lt. assret1(it,1,ip))
			assret(it,1:BLbind,ip)	= agridret(it+1,1,ip)
			conret(it,1:BLbind,ip)	= xgridret(it,1:BLbind,ip) - assret(it,1:BLbind,ip)*annprem(it)
		    IF (QuadraticPref==1) THEN
	    		mucret(it,1:BLbind,ip)    = qprefa*(qprefb-conret(it,1:BLbind,ip))
		    ELSE
        		mucret(it,1:BLbind,ip)    = conret(it,1:BLbind,ip)**(-gam)
		    END IF
		END IF
		
		!interpolate muc1 as fun of ass1 to get muc where BL does not bind
		CALL LinInterp (ngpa,assret1(it,:,ip),conret1(it,:,ip),&
						ngpa-BLbind,agridret(it,BLbind+1:ngpa,ip),conret(it,BLbind+1:ngpa,ip))
        IF (QuadraticPref==1) THEN
            mucret(it,BLbind+1:ngpa,ip) = qprefa*(qprefb - conret(it,BLbind+1:ngpa,ip))
		ELSE
            mucret(it,BLbind+1:ngpa,ip) = conret(it,BLbind+1:ngpa,ip)**(-gam)		
		END IF
		assret(it,BLbind+1:ngpa,ip)	= (xgridret(it,BLbind+1:ngpa,ip) - conret(it,BLbind+1:ngpa,ip))/annprem(it)

		if (any(isnan(conret(it,:,ip)))) THEN
			write(*,*) 'nan encountered'
		end if

	END DO
	!$OMP END PARALLEL DO 
	
END DO	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
!! Working Agents: Rules depend on shocks  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO it = Twork,1,-1
	if (Display==1) write(*,*) 'Solving for decision rules at age ', it				

	!$OMP PARALLEL DO PRIVATE(ie,ip,emuc,iz2,ie2,BLbind,im,lnextm,lpmnext,imnext)
    DO iz = 1,ngpz
	DO ie = 1,ngpe
	DO im = 1,ngpm
	
	    ip = pind(im,iz,ie)
	    
        !solve on tt+1 grid
		IF (it==Twork) THEN			
			muc1(it,:,im,iz,ie) = bet*Rnet* mucret(1,:,pind(im,iz,ie))		!euler equation
        ELSE
			emuc(:) = 0.0
			DO iz2 = 1,ngpz
			DO ie2 = 1,ngpe
		        !find two probabilities and average over:
		        lnextm = (it*mgrid(it,im) + min(ypregrid(it+1,iz2,ie2),pencap))/real(it+1)
			    CALL FindLinProb1 (lnextm,mgrid(it+1,:),imnext,lpmnext)
			    
			    emuc = emuc + lpmnext(1)*muc(it+1,:,imnext(1),iz2,ie2)*ztrans(it,iz,iz2)*edist(it+1,ie2) &
			                + lpmnext(2)*muc(it+1,:,imnext(2),iz2,ie2)*ztrans(it,iz,iz2)*edist(it+1,ie2)
			                
         	END DO
			END DO
			muc1(it,:,im,iz,ie)	=  bet*Rnet*emuc 				!euler equation
		END IF
        
        IF (QuadraticPref==1) THEN
            con1(it,:,im,iz,ie) = qprefb - muc1(it,:,im,iz,ie)/qprefa
        ELSE
            con1(it,:,im,iz,ie) = muc1(it,:,im,iz,ie)**(-1.0/gam) 
        END IF
        
		if (any(isnan(con1(it,:,im,iz,ie)))) THEN
			write(*,*) 'nan encountered'
		end if
						
        IF (it==Twork) THEN	
            ass1(it,:,im,iz,ie) = (con1(it,:,im,iz,ie) +agridret(1,:,ip) - ygrid(it,iz,ie) - tran(it,:,iz,ie))/Rnet  
        ELSE
            ass1(it,:,im,iz,ie) = (con1(it,:,im,iz,ie) +agrid(it+1,:) - ygrid(it,iz,ie) - tran(it,:,iz,ie))/Rnet
        END IF
        
        !deal with borrowing limits
        IF (MINVAL(agrid(it,:)) >= ass1(it,1,im,iz,ie)) THEN		!BL does not bind anywhere
	        BLbind = 0	
        ELSE
            !find point in tt grid where BL starts to bind
	        BLbind = MAXLOC(agrid(it,:), DIM=1,MASK = agrid(it,:) .lt. ass1(it,1,im,iz,ie))						
		    IF (it==Twork) THEN	
			    ass(it,1:BLbind,im,iz,ie) = agridret(1,1,ip)
		    ELSE
			    ass(it,1:BLbind,im,iz,ie) = agrid(it+1,1)
		    END IF
            con(it,1:BLbind,im,iz,ie) = xgrid(it,1:BLbind,iz,ie)-ass(it,1:BLbind,im,iz,ie)
            
            IF (QuadraticPref==1) THEN
                muc(it,1:BLbind,im,iz,ie) = qprefa*(qprefb - con(it,1:BLbind,im,iz,ie) )            
            ELSE
                muc(it,1:BLbind,im,iz,ie) = con(it,1:BLbind,im,iz,ie)**(-gam)
            END IF
        END IF

    					
	    !interpolate muc1 as fun of ass1 to get muc where BL does not bind		
	    CALL LinInterp (ngpa,ass1(it,:,im,iz,ie),con1(it,:,im,iz,ie),&
					    ngpa-BLbind,agrid(it,BLbind+1:ngpa),con(it,BLbind+1:ngpa,im,iz,ie) )
        IF (QuadraticPref==1) THEN
            muc(it,BLbind+1:ngpa,im,iz,ie) = qprefa*(qprefb - con(it,BLbind+1:ngpa,im,iz,ie))
        ELSE
            muc(it,BLbind+1:ngpa,im,iz,ie) = con(it,BLbind+1:ngpa,im,iz,ie)**(-gam)        
        END IF
    
        ass(it,BLbind+1:ngpa,im,iz,ie) = xgrid(it,BLbind+1:ngpa,iz,ie) -con(it,BLbind+1:ngpa,im,iz,ie)	!budget constraint

    END DO  !mgrid
    END DO	!egrid
	END DO	!zgrid
    !$OMP END PARALLEL DO



END DO	!t



END SUBROUTINE Decisions
