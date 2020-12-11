!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine clustering_whole(xl,yl,zl,n,ntime,ncluster,xclust,yclust,zclust,fclust,rms, &
       rmsclust,zrms)
  !                      i  i  i  i   o      o      o      o     o
  !   o      o
  !*****************************************************************************
  !                                                                            *
  !   This routine clusters the particle position into ncluster custers.       *
  !   Input are the longitudes (xl) and latitudes (yl) of the individual       *
  !   points, output are the cluster mean positions (xclust,yclust).           *
  !   Vertical positions are not directly used for the clustering.             *
  !                                                                            *
  !   For clustering, the procedure described in Dorling et al. (1992) is used.*
  !                                                                            *
  !   Dorling, S.R., Davies, T.D. and Pierce, C.E. (1992):                     *
  !   Cluster analysis: a technique for estimating the synoptic meteorological *
  !   controls on air and precipitation chemistry - method and applications.   *
  !   Atmospheric Environment 26A, 2575-2581.                                  *
  !                                                                            *
  !                                                                            *
  !     Author: A. Stohl                                                       *
  !                                                                            *
  !     1 February 2002                                                        *
  !                                                                            *
  ! Variables:                                                                 *
  ! fclust          fraction of particles belonging to each cluster            *
  ! ncluster        number of clusters to be used                              *
  ! rms             total horizontal rms distance after clustering             *
  ! rmsclust        horizontal rms distance for each individual cluster        *
  ! zrms            total vertical rms distance after clustering               *
  ! xclust,yclust,  Cluster centroid positions                                 *
  ! zclust                                                                     *
  ! xl,yl,zl        particle positions                                         *
  !                                                                            *
  !*****************************************************************************
  ! In the R Code  the input ot the distance calculation is a matrix which has one 
  ! row for time along the trajectory and one column for each trajectory.  
  implicit none
  
  

  integer,parameter :: maxpart=100000
  real,parameter :: pi=3.14159265, r_earth=6.371e6, r_air=287.05, ga=9.81
  real,parameter :: cpa=1004.6, kappa=0.286, pi180=pi/180., vonkarman=0.4

  integer :: ntime
  integer :: n,i,j,l,t,nclust(maxpart),ncl,ncluster
  integer :: numb(ncluster)
  real :: xl(ntime,n),yl(ntime,n),zl(ntime,n),xclust(ntime,ncluster),yclust(ntime,ncluster),x(ntime),y(ntime),z(ntime)
  real :: zclust(ntime,ncluster),distance2,distances,distancemin,rms,rmsold
  real :: xav(ntime,ncluster),yav(ntime,ncluster),zav(ntime,ncluster),fclust(ncluster)
  real :: rmsclust(ncluster)
  real :: zdist,zrms
    
  if (n.lt.ncluster) return
  rmsold=-5.

  ! Convert longitude and latitude from degrees to radians
  !*******************************************************

  do i=1,n
    do t=1,ntime
      nclust(i)=i
      xl(t,i)=xl(t,i)*pi180
      yl(t,i)=yl(t,i)*pi180
    end do
  end do


  ! Generate a seed for each cluster
  !*********************************

  do j=1,ncluster
    do t=1,ntime
      zclust(t,j)=0.
      xclust(t,j)=xl(t,j*n/ncluster)
      yclust(t,j)=yl(t,j*n/ncluster)
    end do
  end do


  ! Iterative loop to compute the cluster means
  !********************************************

  do l=1,100

  ! Assign each particle to a cluster: criterion minimum distance to the
  ! cluster mean position
  !*********************************************************************


    do i=1,n
      distancemin=10.**10.
      do j=1,ncluster
        distances=0
        do t=1,ntime 
          distances=distances+distance2(yl(t,i),xl(t,i),yclust(t,j),xclust(t,j))
          if (distances.lt.distancemin) then
            distancemin=distances
            ncl=j 
          endif
        end do
      end do
      nclust(i)=ncl
    end do


  ! Recalculate the cluster centroid position: convert to 3D Cartesian coordinates,
  ! calculate mean position, and re-project this point onto the Earth's surface
  !*****************************************************************************

    do j=1,ncluster
      do t=1, ntime
        xav(t,j)=0.
        yav(t,j)=0.
        zav(t,j)=0.
        rmsclust(j)=0.
        numb(j)=0
      end do
    end do
    rms=0.
    

    do i=1,n
      distances=0.
      numb(nclust(i))=numb(nclust(i))+1
      do t=1,ntime
        distances=distances+distance2(yl(t,i),xl(t,i), &
           yclust(ntime,nclust(i)),xclust(t,nclust(i)))
      end do
  ! rms is the total rms of all particles
  ! rmsclust is the rms for a particular cluster
  !*********************************************

      rms=rms+distances*distances
      rmsclust(nclust(i))=rmsclust(nclust(i))+distances*distances

  ! Calculate Cartesian 3D coordinates from longitude and latitude
  !***************************************************************
      do t=1,ntime

        x(t) = cos(yl(t,i))*sin(xl(t,i))
        y(t) = -1.*cos(yl(t,i))*cos(xl(t,i))
        z(t) = sin(yl(t,i))
        xav(t,nclust(i))=xav(t,nclust(i))+x(t)
        yav(t,nclust(i))=yav(t,nclust(i))+y(t)
        zav(t,nclust(i))=zav(t,nclust(i))+z(t)
      end do
    end do
    
    rms=sqrt(rms/real(n))


  ! Find the mean location in Cartesian coordinates
  !************************************************

    do j=1,ncluster
      do t=1, ntime
        if (numb(j).gt.0) then
          rmsclust(j)=sqrt(rmsclust(j)/real(numb(j)))
          xav(t,j)=xav(t,j)/real(numb(j))
          yav(t,j)=yav(t,j)/real(numb(j))
          zav(t,j)=zav(t,j)/real(numb(j))
        

    ! Project the point back onto Earth's surface
    !********************************************

          xclust(t,j)=atan2(xav(t,j),-1.*yav(t,j))
          yclust(t,j)=atan2(zav(t,j),sqrt(xav(t,j)*xav(t,j)+yav(t,j)*yav(t,j)))
        endif
      end do
    end do


  ! Leave the loop if the RMS distance decreases only slightly between 2 iterations
  !*****************************************************************************

    if ((l.gt.1).and.(abs(rms-rmsold)/rmsold.lt.0.005)) goto 99
    rmsold=rms

  end do

99   continue

  ! Convert longitude and latitude from radians to degrees
  !*******************************************************

  do i=1,n
    do t=1,ntime
      xl(t,i)=xl(t,i)/pi180
      yl(t,i)=yl(t,i)/pi180
      zclust(t,nclust(i))=zclust(t,nclust(i))+zl(t,i)
    end do
  end do

  do j=1,ncluster
    do t=1,ntime
      xclust(t,j)=xclust(t,j)/pi180
      yclust(t,j)=yclust(t,j)/pi180
      if (numb(j).gt.0) zclust(t,j)=zclust(t,j)/real(numb(j))
        fclust(j)=100.*real(numb(j))/real(n)
    end do
  end do

  ! Determine total vertical RMS deviation
  !***************************************

  zrms=0.
  do i=1,n
    zdist=0.
    do t=1,ntime
      zdist=zdist+(zl(t,i)-zclust(t,nclust(i)))
    end do
    zrms=zrms+zdist*zdist
  end do
  if (zrms.gt.0.) zrms=sqrt(zrms/real(n))

end subroutine clustering_whole


function distance2(rlat1,rlon1,rlat2,rlon2)

  !$$$  SUBPROGRAM DOCUMENTATION BLOCK
  !
  ! SUBPROGRAM:  GCDIST     COMPUTE GREAT CIRCLE DISTANCE
  !   PRGMMR: IREDELL       ORG: W/NMC23       DATE: 96-04-10
  !
  ! ABSTRACT: THIS SUBPROGRAM COMPUTES GREAT CIRCLE DISTANCE
  !      BETWEEN TWO POINTS ON THE EARTH. COORDINATES ARE GIVEN IN RADIANS!
  !
  ! PROGRAM HISTORY LOG:
  !   96-04-10  IREDELL
  !
  ! USAGE:    ...GCDIST(RLAT1,RLON1,RLAT2,RLON2)
  !
  !   INPUT ARGUMENT LIST:
  !rlat1    - REAL LATITUDE OF POINT 1 IN RADIANS
  !rlon1    - REAL LONGITUDE OF POINT 1 IN RADIANS
  !rlat2    - REAL LATITUDE OF POINT 2 IN RADIANS
  !rlon2    - REAL LONGITUDE OF POINT 2 IN RADIANS
  !
  !   OUTPUT ARGUMENT LIST:
  !distance2   - REAL GREAT CIRCLE DISTANCE IN KM
  !
  ! ATTRIBUTES:
  !   LANGUAGE: Fortran 90
  !
  !$$$

  use par_mod, only: dp

  implicit none

  real                    :: rlat1,rlon1,rlat2,rlon2,distance2
  real(kind=dp)           :: clat1,clat2,slat1,slat2,cdlon,crd
  real(kind=dp),parameter :: rerth=6.3712e6_dp
  real(kind=dp),parameter :: pi=3.14159265358979_dp

  if ((abs(rlat1-rlat2).lt.0.0003).and. &
       (abs(rlon1-rlon2).lt.0.0003)) then
    distance2=0.0_dp
  else

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    clat1=cos(real(rlat1,kind=dp))
    slat1=sin(real(rlat1,kind=dp))
    clat2=cos(real(rlat2,kind=dp))
    slat2=sin(real(rlat2,kind=dp))
    cdlon=cos(real(rlon1-rlon2,kind=dp))
    crd=slat1*slat2+clat1*clat2*cdlon
    distance2=real(rerth*acos(crd)/1000.0_dp)
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function distance2