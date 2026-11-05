c-----------------------------------------------------------------------
      subroutine cdscal (igeom)
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      LOGICAL          IFCONV

      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT)
     $              ,TB(LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $              ,H2(LX1,LY1,LZ1,LELT)

      include 'ORTHOT'

      if (ifdgfld(ifield)) then
         call cdscal_dg(igeom)
         return
      endif


      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt

      nel    = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         call makeq
         call lagscal

      else                   ! geometry at t^n

         IF (IFPRINT) THEN
            IF (IFIELD.EQ.2.AND.NID.EQ.0)
     $          WRITE (6,*) ' Temperature/Passive scalar solution'
         ENDIF

         if1=ifield-1
         write(name4t,1) if1-1
    1    format('PS',i2)
         if(ifield.eq.2) write(name4t,'(A4)') 'TEMP'

C
C        New geometry
C
         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)
         call bcneusc (ta,-1)
         call add2    (h2,ta,n)
         call add2    (h2,adq(1,1,1,1,ifield-1),n)
         call bcdirsc (t(1,1,1,1,ifield-1))
         call axhelm  (ta,t(1,1,1,1,ifield-1),h1,h2,imesh,ISD)
         call sub3    (tb,bq(1,1,1,1,ifield-1),ta,n)
         call bcneusc (ta,1)
         call add2    (tb,ta,n)

c        PHYSALIS DIAGNOSTIC (added by explicit user instruction; original
c        file preserved unmodified as conduct_ori.f alongside this one).
c        tb here is the REAL solver RHS residual (bq-A*t, post-bcneusc),
c        captured at the exact instant core/conduct.f itself computes it,
c        before hsolve/masking -- this is the quantity a .usr-level
c        post-hoc reconstruction of axhelm(t)-bq could not reproduce
c        (confirmed nonzero ~-97 "free fluid" residual in that attempt).
         if (ifaxis.and.ifaziv.and.ifield.eq.2)
     $      call physalis_swirl_reaction_diag(tb,iter,nel)

         if(iftmsh(ifield)) then
           call hsolve  (name4t,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),bintm1)
         else
           call hsolve  (name4t,TA,TB,H1,H2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         endif 

         call add2    (t(1,1,1,1,ifield-1),ta,n)
         call physalis_restore_swirl_cage_tmask(ifield)
         call physalis_sync_axisym_swirl(ifield)

c        PHYSALIS DIAGNOSTIC, POST-SOLVE (added by explicit user
c        instruction; original file preserved as conduct_ori.f).
c        The pre-solve tb split (physalis_swirl_reaction_diag above)
c        is a load/RHS split, not a discrete reaction -- this instead
c        recomputes the residual axhelm(t)-bq using the FINAL, post-
c        update/post-mask-restore/post-sync t field and the SAME h1,h2
c        this iteration's solve used, so a nonzero value at cage nodes
c        is the actual discrete Dirichlet-constraint reaction and a
c        near-zero value at free nodes confirms the solve is truly
c        converged there.
         if (ifaxis.and.ifaziv.and.ifield.eq.2)
     $      call physalis_swirl_reaction_post_diag(h1,h2,iter,nel,ISD)

         call cvgnlps (ifconv) ! Check convergence for nonlinear problem
         if (ifconv) goto 2000

C        Radiation case, smooth convergence, avoid flip-flop (ER).
         call cmult (ta,0.5,n)
         call sub2  (t(1,1,1,1,ifield-1),ta,n)

 1000    continue
 2000    continue

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine makeuq

c     Fill up user defined forcing function and collocate will the
c     mass matrix on the Gauss-Lobatto mesh.

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'

      n = lx1*ly1*lz1*nelfld(ifield)

      if (.not.ifcvfld(ifield)) time = time-dt ! Set time to t^n-1 for user function

      if (nio.eq.0.and.loglevel.gt.2) 
     $   write(6,*) 'makeuq', ifield, time
      call setqvol(bq(1,1,1,1,ifield-1),adq(1,1,1,1,ifield-1))
      call col2   (bq(1,1,1,1,ifield-1) ,bm1,n)

      if (.not.ifcvfld(ifield)) time = time+dt ! Restore time

      return
      end
c-----------------------------------------------------------------------
      subroutine setqvol(bql,aql)

c     Set user specified volumetric forcing function (e.g. heat source).

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'

      real bql(lx1*ly1*lz1,lelt),aql(lx1*ly1*lz1,lelt)

#ifdef TIMER
      etime1=dnekclock()
#endif

      nel   = nelfld(ifield)
      nxyz1 = lx1*ly1*lz1
      n     = nxyz1*nel

      do iel=1,nel

         call nekuq (bql,aql,iel) ! ONLY SUPPORT USERQ - pff, 3/08/16

c        igrp = igroup(iel)
c        if (matype(igrp,ifield).eq.1) then ! constant source within a group
c           cqvol = cpgrp(igrp,ifield,3)
c           call cfill (bql(1,iel),cqvol,nxyz1)
c        else  !  pff 2/6/96 ............ default is to look at userq
c           call nekuq (bql,iel)
c        endif

      enddo
c
c 101 FORMAT(' Wrong material type (',I3,') for group',I3,', field',I2
c    $    ,/,' Aborting in SETQVOL.')
C  

#ifdef TIMER
      tusfq=tusfq+(dnekclock()-etime1)
#endif
 
      return
      end
C
      subroutine nekuq (bql,aql,iel)
C------------------------------------------------------------------
C
C     Generate user-specified volumetric source term (temp./p.s.)
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'PARALLEL'
      include 'TSTEP'
      include 'NEKUSE'
      include 'INPUT'
c
      real bql(lx1,ly1,lz1,lelt),aql(lx1,ly1,lz1,lelt)
c
      ielg = lglel(iel)
      do 10 k=1,lz1
      do 10 j=1,ly1
      do 10 i=1,lx1
         if (optlevel.le.2) call nekasgn (i,j,k,iel)
         qvol = 0.0
         avol = 0.0
         call userq   (i,j,k,ielg)
         bql(i,j,k,iel) = qvol
         aql(i,j,k,iel) = avol
 10   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine convab
C---------------------------------------------------------------
C
C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C
C---------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      common /scruz/ ta (lx1*ly1*lz1*lelt)

      nel = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      call convop  (ta,t(1,1,1,1,ifield-1))
      do i=1,n
        bq(i,1,1,1,ifield-1) = bq (i,1,1,1,ifield-1)
     $                       - bm1(i,1,1,1)*ta(i)*vtrans(i,1,1,1,ifield)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makeabq
C
C     Sum up contributions to 3rd order Adams-Bashforth scheme.
C
      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'

      ab0   = ab(1)
      ab1   = ab(2)
      ab2   = ab(3)
      nel   = nelfld(ifield)
      n     = lx1*ly1*lz1*nel

      do i=1,n
         ta=ab1*vgradt1(i,1,1,1,ifield-1)+ab2*vgradt2(i,1,1,1,ifield-1)
         vgradt2(i,1,1,1,ifield-1)=vgradt1(i,1,1,1,ifield-1)
         vgradt1(i,1,1,1,ifield-1)=bq     (i,1,1,1,ifield-1)
         bq     (i,1,1,1,ifield-1)=bq     (i,1,1,1,ifield-1)*ab0+ta
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine makebdq
C-----------------------------------------------------------------------
C
C     Add contributions to F from lagged BD terms.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ tb(lt),h2(lt)

      nel   = nelfld(ifield)
      n     = lx1*ly1*lz1*nel
   
      const = 1./dt
      do i=1,n
         h2(i)=const*vtrans(i,1,1,1,ifield)
         tb(i)=bd(2)*bm1(i,1,1,1)*t(i,1,1,1,ifield-1)
      enddo

      do ilag=2,nbd
         if (ifgeom) then
            do i=1,n
               ta=bm1lag(i,1,1,1,ilag-1)*tlag(i,1,1,1,ilag-1,ifield-1)
               tb(i)=tb(i)+ta*bd(ilag+1)
            enddo
         else
            do i=1,n
               ta=bm1(i,1,1,1)*tlag(i,1,1,1,ilag-1,ifield-1)
               tb(i)=tb(i)+ta*bd(ilag+1)
            enddo
         endif
      enddo

      call addcol3 (bq(1,1,1,1,ifield-1),tb,h2,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine makebdq_solid
C-----------------------------------------------------------------------
C
C     Add contributions to F from lagged BD terms.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ tb(lt),h2(lt)

      nel   = nelfld(ifield)
      n     = lx1*ly1*lz1*nelv
      ndif=lx1*ly1*lz1*(nel-nelv)

      if(ndif.eq.0) return

      const = 1./dt
      do i=1,ndif
         h2(n+i)=const*vtrans(n+i,1,1,1,ifield)
         tb(n+i)=bd(2)*bm1(n+i,1,1,1)*t(n+i,1,1,1,ifield-1)
      enddo

      do ilag=2,nbd
         if (ifgeom) then
            do i=1,ndif
               ta=bm1lag(n+i,1,1,1,ilag-1)
     $             *tlag(n+i,1,1,1,ilag-1,ifield-1)
               tb(n+i)=tb(n+i)+ta*bd(ilag+1)
            enddo
         else
            do i=1,ndif
               ta=bm1(n+i,1,1,1)*tlag(n+i,1,1,1,ilag-1,ifield-1)
               tb(n+i)=tb(n+i)+ta*bd(ilag+1)
            enddo
         endif
      enddo

      call addcol3 (bq(n+1,1,1,1,ifield-1),tb(n+1),h2(n+1),ndif)

      return
      end
c-----------------------------------------------------------------------
      subroutine lagscal   !  Keep old passive scalar field(s) 

      include 'SIZE'
      include 'TOTAL'

      n = lx1*ly1*lz1*nelfld(ifield)

      do ilag=nbdinp-1,2,-1
         call copy (tlag(1,1,1,1,ilag  ,ifield-1),
     $              tlag(1,1,1,1,ilag-1,ifield-1),n)
      enddo

      call copy (tlag(1,1,1,1,1,ifield-1),t(1,1,1,1,ifield-1),n)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfldrq (x,txt10,ichk)
      include 'SIZE'
      include 'TSTEP'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /-3/

      if (idum.lt.0) return
c
C
      mtot = lx1*ly1*lz1*nelv
      if (lx1.gt.8.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      nell = nelt
      rnel = nell
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nell-ne+1
      k = 1
      do ie=1,1
         ne = 0
         write(6,116) txt10,k,ie,xmin,xmax,istep,time
         do l=0,1
            write(6,117) 
            do j=ly1,1,-1
              if (lx1.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
              if (lx1.eq.8) write(6,118) ((x(i,j,k,e+l),i=1,lx1),e=1,1)
            enddo
         enddo
      enddo
C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
  118 FORMAT(8f12.9)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine cdscal_expl (igeom)
C
C     explicit convection-diffusion equation for passive scalar
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      common  /cprint/ ifprint
      logical          ifprint
      logical          ifconv

      common /scrns/ ta(lx1,ly1,lz1,lelt)
     $              ,tb(lx1,ly1,lz1,lelt)
      common /scrvh/ h1(lx1,ly1,lz1,lelt)
     $              ,h2(lx1,ly1,lz1,lelt)


c     QUESTIONABLE support for Robin BC's at this point! (5/15/08)

      nel    = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}

         call makeq
         call lagscal

      else                   ! geometry at t^n

         if (nio.eq.0) write (6,*) istep,ifield,' explicit step'


C        New geometry

         isd = 1
         if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)

         call bcneusc (ta,-1)       ! Modify diagonal for Robin condition
         call add2    (h2,ta ,n)
         call col2    (h2,BM1,n)

         call bcneusc (tb,1)        ! Modify rhs for flux bc
         call add2    (bq(1,1,1,1,ifield-1),tb,n)

         call dssum   (bq(1,1,1,1,ifield-1),lx1,ly1,lz1)
         call dssum   (h2,lx1,ly1,lz1)

         call invcol3 (t(1,1,1,1,ifield-1),bq(1,1,1,1,ifield-1),h2,n)

         call bcdirsc (t(1,1,1,1,ifield-1)) ! --> no mask needed

      endif                   ! geometry at t^n

      return
      end
c-----------------------------------------------------------------------
      subroutine diffab  ! explicit treatment of diffusion operator
c
c     Eulerian scheme, add diffusion term to forcing function 
c     at current time step.
c

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'INPUT'

      common /scruz/ ta(lx1,ly1,lz1,lelt)
     $              ,h2(lx1,ly1,lz1,lelt)

      nel = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      intype = 0
      if (iftran) intype = -1

      isd = 1
      if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

      imesh = 1
c      if (iftmsh(ifield)) imesh=2

      call rzero   (h2,n)
      call axhelm  (ta,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield)
     $             ,h2,imesh,isd)
      call sub2    (bq(1,1,1,1,ifield-1),ta,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_eta_alpha2

c     Set up required dg terms, e.g., alpha, eta, etc.
c     Face weight: .5 interior, 1. boundary

      include 'SIZE'
      include 'TOTAL'
      common /mysedg/ eta


      integer e,f,pf

      nface = 2*ldim
      call dsset(lx1,ly1,lz1)

      eta =  5         !   Semi-optimized value, single domain

      do e=1,nelfld(ifield)
      do f=1,nface

         pf     = eface1(f)
         js1    = skpdat(1,pf)
         jf1    = skpdat(2,pf)
         jskip1 = skpdat(3,pf)
         js2    = skpdat(4,pf)
         jf2    = skpdat(5,pf)
         jskip2 = skpdat(6,pf)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1
            i = i+1
            a = area(i,1,f,e)  ! Check Fydkowski notes 
            a = a*a*fw(f,e)    ! For ds_avg used below, plus quad weight
            etalph(i,f,e) = eta*(a/bm1(j1,j2,1,e))
c           write(6,*) i,j1,j2,e,f,a,etalph(i,f,e)
         enddo
         enddo
      enddo
      enddo

      call fgslib_gs_op (dg_hndlx,etalph,1,1,0)  ! 1 ==> +

      return
      end
c-----------------------------------------------------------------------
      subroutine fwght(msk1,mult)
      include 'SIZE'
      include 'TOTAL'
      parameter (lx=lx1*ly1*lz1)
      real msk1(lx1,ly1,lz1),mult(lx1,ly1,lz1,1)
      integer e,f

      parameter(lf=lx1*lz1*2*ldim*lelt)
      common /scrdg/uf(lx1*lz1,2*ldim,lelt)



      do i=1,lf
         uf(i,1,1)=1.
      enddo
      call fgslib_gs_op (dg_hndlx,uf,1,1,0)  ! 1 ==> +

      nface = 2*ldim
      do e=1,nelt
      do f=1,nface
         fw(f,e) = 1                        ! Boundary
         if (uf(1,f,e).gt.1.1) fw(f,e)=0.5  ! Interior
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_setup
      include 'SIZE'
      include 'TOTAL'

      common /ivrtx/ vertex ((2**ldim)*lelt)
      common /ctmp1/ qs(lx1*ly1*lz1*lelt)
      integer*8 vertex

      logical ifany
      save    ifany
      data    ifany /.false./

      do ifield=1,ldimt1
         if (ifdgfld(ifield)) ifany=.true.
      enddo

      if (ifany) then

        call setup_dg_gs(dg_hndlx,lx1,ly1,lz1,nelt,nelgt,vertex)
        call dg_set_fc_ptr

        param(59)=1
        call geom_reset(1)
        call set_unr
      endif


      n = lx1*ly1*lz1*nelt
      call invers2(binvdg,bm1,n)

      call set_dg_wgts

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_setup2(mask)
      include 'SIZE'
      include 'TOTAL'

      real mask(1)
      common /ctmp0/ qs(lx1*ly1*lz1*lelt)

      integer ifld_last
      save    ifld_last
      data    ifld_last /0/

      if (ifield.eq.ifld_last) return
      ifld_last = ifield

      call fwght  (qs,mask)  ! qs is work array
      call set_eta_alpha2     ! Set eta/h

      return
      end
c-----------------------------------------------------------------------
      subroutine cdscal_dg (igeom)
C
C     Solve the convection-diffusion equation for passive scalar IPSCAL
C
      include 'SIZE'
      include 'TOTAL'
      common  /cprint/ ifprint
      logical          ifprint,ifconv

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ ta(lt),tb(lt)
      common /scrvh/ h1(lt),h2(lt)

      include 'ORTHOT'  ! This must be fixed

      call dg_setup2(tmask(1,1,1,1,ifield-1))

      nel = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      if (igeom.eq.1) then   ! old geometry at t^{n-1}

         call makeq
         call lagscal

      else                   ! new geometry at t^n

         if (ifprint.and.nio.eq.0) 
     $      write (6,*) ' Temperature/Passive scalar solution',ifield

         if1=ifield-1
         write(name4t,1) if1-1
    1    format('PS',i2)
         if (ifield.eq.2) write(name4t,'(A4)') 'TEMP'
 
         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)

c        call bcneusc (ta,-1)     !! Not Yet supported for DG
c        call add2    (h2,ta,n)   !! Not Yet supported for DG

         call rzero             (tb,n)
         call bcdirsc           (   t (1,1,1,1,ifield-1))
         call conv_bdry_dg_weak (tb,t (1,1,1,1,ifield-1))
         call hxdg_surfa        (tb,t (1,1,1,1,ifield-1),h1,h2)
         call add2              (tb,bq(1,1,1,1,ifield-1),n)

         call hmholtz_dg(name4t,t(1,1,1,1,ifield-1),tb,h1,h2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tolht(ifield),nmxt(ifield-1))
         call physalis_restore_swirl_cage_tmask(ifield)
         call physalis_sync_axisym_swirl(ifield)

      endif  ! End of IGEOM branch.

      return
      end
c-----------------------------------------------------------------------
      subroutine physalis_swirl_reaction_diag(tb,iter,nel)
c  PHYSALIS DIAGNOSTIC (added by explicit user instruction; original
c  file preserved unmodified as conduct_ori.f alongside this one).
c  Sums the REAL swirl-solve RHS residual tb (=bq-A*t, post-bcneusc,
c  pre-hsolve/pre-mask -- see the call site in cdscal above) times
c  nodal radius, separately over (a) the velocity-cage Dirichlet-
c  constrained swirl nodes (mask_cage_u==1) and (b) ordinary free-
c  fluid nodes (neither cage nor inside-particle). tb is already a
c  quadrature-weighted weak-form load (same units as bq), so summing
c  R*tb directly gives a torque-like quantity with no extra mass-
c  matrix weighting needed.
c
c  T_cage_tb/T_free_tb: as computed (tb sign convention, =bq-A*t).
c  T_cage_minus_tb/T_free_minus_tb: the negation, reported alongside
c  since the correct sign for "reaction torque supplied by the
c  solver to hold the Dirichlet value" was not assumed in advance --
c  compare both against gap_Ld_minus_Lp from the .usr-side
c  physalis_axisym_angular_budget.dat to see which sign matches.
c  T_free_tb (and its negation) is a SANITY CHECK: this should be
c  close to zero at a converged nonlinear iteration, since ordinary
c  free nodes have no Dirichlet constraint forcing tb away from what
c  hsolve would otherwise drive to zero -- unlike the earlier .usr-
c  level post-hoc reconstruction attempt, which gave T_interior_
c  check_z~=-97 (not small), because it could not reproduce the exact
c  h1/h2/bq/mask state at the true solve instant. This version reads
c  tb at that exact instant, inside cdscal itself.
      include 'SIZE'
      include 'TOTAL'

      real tb(lx1,ly1,lz1,lelt)
      integer iter,nel

      integer mask_inside_u(lx1,ly1,lz1,lelv)
      integer mask_cage_u(lx1,ly1,lz1,lelv)
      common /physalis_masks_u/ mask_inside_u,mask_cage_u

      integer sp_cons_mode
      common /sp_cons_i/ sp_cons_mode

      real Tcage_loc,Tfree_loc,rabs
      real Tcage_g,Tfree_g
      integer n_cage_loc,n_free_loc,n_cage_g,n_free_g
      integer e,i,j,k,fu,write_stride
      parameter(fu=97,write_stride=100)
      logical first_rd
      save first_rd
      data first_rd/.true./

      if (sp_cons_mode.ne.2) return
      if (istep.gt.10 .and. mod(istep,write_stride).ne.0) return

      Tcage_loc = 0.0
      Tfree_loc = 0.0
      n_cage_loc = 0
      n_free_loc = 0
      do e=1,nelv
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         rabs = ym1(i,j,k,e)
         if (mask_cage_u(i,j,k,e).eq.1) then
            Tcage_loc = Tcage_loc + rabs*tb(i,j,k,e)
            n_cage_loc = n_cage_loc + 1
         else if (mask_inside_u(i,j,k,e).ne.1) then
            Tfree_loc = Tfree_loc + rabs*tb(i,j,k,e)
            n_free_loc = n_free_loc + 1
         endif
      enddo
      enddo
      enddo
      enddo

      Tcage_g = glsum(Tcage_loc,1)
      Tfree_g = glsum(Tfree_loc,1)
      n_cage_g = iglsum(n_cage_loc,1)
      n_free_g = iglsum(n_free_loc,1)

      if (nid.eq.0) then
         if (first_rd) then
            open(fu,file='physalis_axisym_swirl_reaction_core.dat',
     $           status='replace')
            write(fu,'(A)') '# step time iter T_cage_tb '//
     $       'T_cage_minus_tb T_free_tb T_free_minus_tb '//
     $       'n_cage n_free'
            first_rd = .false.
         else
            open(fu,file='physalis_axisym_swirl_reaction_core.dat',
     $           status='old',position='append')
         endif
         write(fu,'(I8,1X,E18.10,1X,I4,1X,4E18.10,1X,2I8)')
     $    istep,time,iter,Tcage_g,-Tcage_g,Tfree_g,-Tfree_g,
     $    n_cage_g,n_free_g
         close(fu)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine physalis_swirl_reaction_post_diag(h1,h2,iter,nel,isd)
c  PHYSALIS DIAGNOSTIC, POST-SOLVE (added by explicit user instruction;
c  original file preserved unmodified as conduct_ori.f alongside this
c  one).
c
c  physalis_swirl_reaction_diag (above) sums the PRE-solve RHS split
c  tb=bq-A*t -- a load/RHS decomposition, not a discrete reaction (the
c  user's own reading of the first round of results: T_cage_tb came
c  out ~8x larger in magnitude than the Ld_z-Lp_z gap, and T_free_tb
c  was NOT small, so tb itself cannot be read as "the reaction").
c
c  This instead recomputes, using the FINAL t field (after add2/
c  restore_tmask/sync_axisym_swirl have all already run this
c  iteration) and the SAME h1,h2 this iteration's hsolve used:
c     res = axhelm(t,h1,h2,imesh,isd) - bq
c  At FREE (non-cage, non-inside) nodes the converged solve should
c  make this close to zero (a direct check of solve consistency, not
c  assumed). At CAGE (Dirichlet-masked) nodes, res is the actual
c  discrete reaction load the solver supplies to hold that value --
c  summing R*res there is the physically meaningful reaction torque
c  the earlier pre-solve tb split could not provide.  T_inside
c  (r<a, rigid-body-only nodes, never Lamb-corrected) is reported too
c  for completeness/symmetry, not because it is expected to matter.
c
c  Does not call any .usr-defined subroutine (only core utilities
c  axhelm/glsum/iglsum and the physalis_masks_u/sp_cons_i common
c  blocks, matching physalis_swirl_reaction_diag's own precedent)
c  so this stays link-safe for every case that builds against this
c  shared core, angular or not.
      include 'SIZE'
      include 'TOTAL'

      real h1(lx1,ly1,lz1,lelt), h2(lx1,ly1,lz1,lelt)
      integer iter,nel,isd

      real resw(lx1,ly1,lz1,lelt)
      common /phys_react_work/ resw

      integer mask_inside_u(lx1,ly1,lz1,lelv)
      integer mask_cage_u(lx1,ly1,lz1,lelv)
      common /physalis_masks_u/ mask_inside_u,mask_cage_u

      integer sp_cons_mode
      common /sp_cons_i/ sp_cons_mode

      real Tcage_loc,Tfree_loc,Tinside_loc,rabs
      real Tcage_g,Tfree_g,Tinside_g
      integer n_cage_loc,n_free_loc,n_inside_loc
      integer n_cage_g,n_free_g,n_inside_g
      integer e,i,j,k,fu,write_stride,npt
      parameter(fu=98,write_stride=100)
      logical first_pd
      save first_pd
      data first_pd/.true./

      if (sp_cons_mode.ne.2) return
      if (istep.gt.10 .and. mod(istep,write_stride).ne.0) return

      npt = lx1*ly1*lz1*nelv
      call axhelm(resw,t(1,1,1,1,ifield-1),h1,h2,imesh,isd)
      call sub2(resw,bq(1,1,1,1,ifield-1),npt)

      Tcage_loc = 0.0
      Tfree_loc = 0.0
      Tinside_loc = 0.0
      n_cage_loc = 0
      n_free_loc = 0
      n_inside_loc = 0
      do e=1,nelv
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         rabs = ym1(i,j,k,e)
         if (mask_cage_u(i,j,k,e).eq.1) then
            Tcage_loc = Tcage_loc + rabs*resw(i,j,k,e)
            n_cage_loc = n_cage_loc + 1
         else if (mask_inside_u(i,j,k,e).eq.1) then
            Tinside_loc = Tinside_loc + rabs*resw(i,j,k,e)
            n_inside_loc = n_inside_loc + 1
         else
            Tfree_loc = Tfree_loc + rabs*resw(i,j,k,e)
            n_free_loc = n_free_loc + 1
         endif
      enddo
      enddo
      enddo
      enddo

      Tcage_g = glsum(Tcage_loc,1)
      Tfree_g = glsum(Tfree_loc,1)
      Tinside_g = glsum(Tinside_loc,1)
      n_cage_g = iglsum(n_cage_loc,1)
      n_free_g = iglsum(n_free_loc,1)
      n_inside_g = iglsum(n_inside_loc,1)

      if (nid.eq.0) then
         if (first_pd) then
            open(fu,file='physalis_axisym_swirl_reaction_post.dat',
     $           status='replace')
            write(fu,'(A)') '# step time iter T_cage_res T_free_res '//
     $       'T_inside_res n_cage n_free n_inside'
            first_pd = .false.
         else
            open(fu,file='physalis_axisym_swirl_reaction_post.dat',
     $           status='old',position='append')
         endif
         write(fu,'(I8,1X,E18.10,1X,I4,1X,3E18.10,1X,3I8)')
     $    istep,time,iter,Tcage_g,Tfree_g,Tinside_g,
     $    n_cage_g,n_free_g,n_inside_g
         close(fu)
      endif

      return
      end
c-----------------------------------------------------------------------
