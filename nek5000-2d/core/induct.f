      subroutine induct (igeom)
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'

      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      ifield = ifldmhd

      if (igeom.eq.1) then

         call makebsource_mhd

      else

         call lagbfield
         call lagvel

         call elsasserh(igeom)

         call vol_flow

      endif
      return
      end
      subroutine lagbfield
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      do ilag=nbdinp-1,2,-1
         call opcopy
     $    ( bxlag(1,ilag  ),bylag(1,ilag  ),bzlag(1,ilag  )
     $    , bxlag(1,ilag-1),bylag(1,ilag-1),bzlag(1,ilag-1) )
      enddo
      call opcopy (bxlag,bylag,bzlag,bx,by,bz)
      return
      end
      subroutine makebsource_mhd
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'CTIMER'

      if (icalld.eq.0) tbmhd=0.0
      icalld = icalld+1
      nbmhd  = icalld
      etime1 = dnekclock()
      ifield = 1
                                      call makeuf
      ifield = ifldmhd
                                      call makeufb
      if (ifaxis) then
         do ifield = 2,npscal+1
                                      call makeuq
         enddo
         ifield = ifldmhd
      endif

      if (ifnav.and.(.not.ifchar)) then
         call advab_elsasser_fast
      endif
      if (ifchar) then
         write(6,*) 'No IFCHAR for MHD, yet.'
         call exitt
      endif

      ifield = 1
      if (iftran)                     call makeabf
                                      call makebdf
      ifield = ifldmhd
      if (iftran)                     call makextb
                                      call makebdfb
      tbmhd=tbmhd+(dnekclock()-etime1)
      return
      end
      subroutine makeufb
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      time = time-dt
      call nekuf   (bmx,bmy,bmz)
      call opcolv2 (bmx,bmy,bmz,vtrans(1,1,1,1,ifield),bm1)
      time = time+dt
      return
      end
      subroutine makextb
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv
      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,bbx1,bbx2,ab1,ab2,ntot1)
      call add3s2 (ta2,bby1,bby2,ab1,ab2,ntot1)
      call copy   (bbx2,bbx1,ntot1)
      call copy   (bby2,bby1,ntot1)
      call copy   (bbx1,bmx,ntot1)
      call copy   (bby1,bmy,ntot1)
      call add2s1 (bmx,ta1,ab0,ntot1)
      call add2s1 (bmy,ta2,ab0,ntot1)
      if (ldim.eq.3) then
         call add3s2 (ta3,bbz1,bbz2,ab1,ab2,ntot1)
         call copy   (bbz2,bbz1,ntot1)
         call copy   (bbz1,bmz,ntot1)
         call add2s1 (bmz,ta3,ab0,ntot1)
      endif
      return
      end
      subroutine makebdfb
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      COMMON /SCRNS/ TA1(LX1,LY1,LZ1,LELV)
     $ ,             TA2(LX1,LY1,LZ1,LELV)
     $ ,             TA3(LX1,LY1,LZ1,LELV)
     $ ,             TB1(LX1,LY1,LZ1,LELV)
     $ ,             TB2(LX1,LY1,LZ1,LELV)
     $ ,             TB3(LX1,LY1,LZ1,LELV)
     $ ,             H2 (LX1,LY1,LZ1,LELV)
      NTOT1 = lx1*ly1*lz1*NELV
      CONST = 1./DT
      CALL CMULT2(H2,vtrans(1,1,1,1,ifield),CONST,NTOT1)
      CALL OPCOLV3c (TB1,TB2,TB3,BX,BY,BZ,BM1,bd(2))
      DO 100 ILAG=2,NBD
         IF (IFGEOM) THEN
            CALL OPCOLV3c(TA1,TA2,TA3,BXLAG (1,ILAG-1),
     $                                BYLAG (1,ILAG-1),
     $                                BZLAG (1,ILAG-1),
     $                                BM1LAG(1,1,1,1,ILAG-1),bd(ilag+1))
         ELSE
            CALL OPCOLV3c(TA1,TA2,TA3,BXLAG (1,ILAG-1),
     $                                BYLAG (1,ILAG-1),
     $                                BZLAG (1,ILAG-1),
     $                                BM1                   ,bd(ilag+1))
         ENDIF
         CALL OPADD2  (TB1,TB2,TB3,TA1,TA2,TA3)
 100  CONTINUE
      CALL OPADD2col (BMX,BMY,BMZ,TB1,TB2,TB3,h2)
      return
      end
      subroutine cresvib(resv1,resv2,resv3,h1,h2)
      include 'SIZE'
      include 'TOTAL'
      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
      call bcdirvc (bx,by,bz,b1mask,b2mask,b3mask)
      call extrapp (pm,pmlag)
      call opgradt (resv1,resv2,resv3,pm)
      call opadd2  (resv1,resv2,resv3,bmx,bmy,bmz)
      call ophx    (w1,w2,w3,bx,by,bz,h1,h2)
      call opsub2  (resv1,resv2,resv3,w1,w2,w3)
      return
      end
      subroutine cresvibp(resv1,resv2,resv3,h1,h2)
      include 'SIZE'
      include 'TOTAL'
      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
      call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)
      if (ifstrs)  call bcneutr
      call extrapp (pr,prlag)
      call opgradt (resv1,resv2,resv3,pr)
      call opadd2  (resv1,resv2,resv3,bfx,bfy,bfz)
      call ophx    (w1,w2,w3,vx,vy,vz,h1,h2)
      call opsub2  (resv1,resv2,resv3,w1,w2,w3)
      return
      end
      subroutine incomprn (ux,uy,uz,up)
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      parameter(nset = 1 + lbelv/lelv)
      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
      common /orthbi/ nprv(2)
      logical ifprjp
      integer phys_mode,n_phys_iter,lmax,pressure_mode
      real beta_coeff,phys_tol_coeff,phys_tol_cage
      real r_sample_over_a,cage_thickness_u,cage_thickness_p
      common /physalis_params/ phys_mode,n_phys_iter,lmax,
     $ pressure_mode,beta_coeff,phys_tol_coeff,phys_tol_cage,
     $ r_sample_over_a,cage_thickness_u,cage_thickness_p

      ifprjp=.false.
      istart=param(95)
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.
      if (phys_mode.ge.5) ifprjp=.false.

      if (icalld.eq.0) tpres=0.0
      icalld = icalld+1
      npres  = icalld
      etime1 = dnekclock()

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      intype = 1

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opdiv   (dp,ux,uy,uz)

      bdti = -bd(1)/dt
      call cmult   (dp,bdti,ntot2)

      call add2col2(dp,bm2,usrdiv,ntot2)

      call stage23_pressure_rhs_ortho(dp)

      i = 1 + ifield/ifldmhd
      if (ifprjp)   call setrhsp  (dp,h1,h2,h2inv,pset(1,i),nprv(i))
                    scaledt = dt/bd(1)
                    scaledi = 1./scaledt
                    call cmult(dp,scaledt,ntot2)
      call stage23_pressure_target(dp,up,scaledt)
      call stage23_operator_audit(h1,h2,h2inv,intype)
      call stage23_pressure_rhs_lift(dp,h1,h2,h2inv,
     $ intype)
                    call esolver  (dp,h1,h2,h2inv,intype)
      if (phys_mode.ge.5) call stage23_pressure_restore(dp)
                    call cmult(dp,scaledi,ntot2)
      if (ifprjp)   call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      call add2(up,dp,ntot2)
      if (phys_mode.ge.5) call physalis_enforce_pressure_cage(up)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)

c     Stage35BI: apply exactly the same free-velocity projector F used
c     by CDABDTP, so the solved Schur operator and applied correction
c     are algebraically consistent.
      call physalis_apply_u_free_mask(dv1,dv2,dv3,1)

      dtb  = dt/bd(1)
      call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )

      if (ifmhd)  call chkptol

      tpres=tpres+(dnekclock()-etime1)

      return
      end
      subroutine extrapp(p,plag)
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'

      real  p    (lx2,ly2,lz2,1)
     $     ,plag (lx2,ly2,lz2,1)

      common /cgeom/ igeom

      ntot2 = lx2*ly2*lz2*nelv

      if (nbd.eq.2.and.nbdinp.gt.2.and.igeom.le.2) then
         call copy(plag,p,ntot2)
      elseif (nbd.eq.3.and.igeom.le.2) then

         const = dtlag(1)/dtlag(2)

         do i=1,ntot2
            pnm1          = p   (i,1,1,1)
            pnm2          = plag(i,1,1,1)
            p   (i,1,1,1) = pnm1 + const*(pnm1-pnm2)
            plag(i,1,1,1) = pnm1
         enddo

      elseif (nbd.gt.3) then
         WRITE (6,*) 'Pressure extrapolation cannot be completed'
         WRITE (6,*) 'Try a lower-order temporal scheme'
         call exitt
      endif
      return
      end
      subroutine opzero(ux,uy,uz)
      include 'SIZE'
      include 'TOTAL'
      real ux(1),uy(1),uz(1)
      n = lx1*ly1*lz1*nelfld(ifield)
      call rzero(ux,n)
      call rzero(uy,n)
      if (if3d) call rzero(uz,n)
      return
      end
      subroutine opnorm(unorm,ux,uy,uz,type3)
      include 'SIZE'
      include 'TOTAL'
      character*3 type3
      real ux(1),uy(1),uz(1)
      real un(3),wn(3)
      n = lx1*ly1*lz1*nelfld(ifield)
      if (type3.eq.'L2 ') then
         if (if3d) then
            un(1) = vlsc3(ux,ux,bm1,n)
            un(2) = vlsc3(uy,uy,bm1,n)
            un(3) = vlsc3(uz,uz,bm1,n)
            un(1) = un(1) + un(2) + un(3)
            unorm = glsum(un(1),1)
            if (unorm.gt.0) unorm = sqrt(unorm/volvm1)
         else
            un(1) = vlsc3(ux,ux,bm1,n)
            un(2) = vlsc3(uy,uy,bm1,n)
            un(1) = un(1) + un(2)
            unorm = glsum(un(1),1)
            if (unorm.gt.0) unorm = sqrt(unorm/volvm1)
         endif
      endif
      return
      end
      subroutine lorentz_force (lf,b1,b2,b3,w1,w2)
      include 'SIZE'
      real lf(lx1*ly1*lz1*lelv,ldim)
      real b1(lx1*ly1*lz1*lelv)
      real b2(lx1*ly1*lz1*lelv)
      real b3(lx1*ly1*lz1*lelv)
      call curl(lf,b1,b2,b3,.false.,w1,w2)
      ntot = lx1*ly1*lz1*nelv
      do i=1,ntot
         c1 = lf(i,2)*b3(i) - lf(i,3)*b2(i)
         c2 = lf(i,3)*b1(i) - lf(i,1)*b3(i)
         c3 = lf(i,1)*b2(i) - lf(i,2)*b1(i)
         lf(i,1) = c1
         lf(i,2) = c2
         lf(i,3) = c3
      enddo
      return
      end
      subroutine curl(vort,u,v,w,ifavg,work1,work2)
      include 'SIZE'
      include 'TOTAL'
      logical ifavg
      parameter(lt=lx1*ly1*lz1*lelv)
      real vort(lt,3),work1(1),work2(1),u(1),v(1),w(1)
      ntot  = lx1*ly1*lz1*nelv
      if (if3d) then
           call dudxyz(work1,w,rym1,sym1,tym1,jacm1,1,2)
           call dudxyz(work2,v,rzm1,szm1,tzm1,jacm1,1,3)
           call sub3(vort(1,1),work1,work2,ntot)
           call dudxyz(work1,u,rzm1,szm1,tzm1,jacm1,1,3)
           call dudxyz(work2,w,rxm1,sxm1,txm1,jacm1,1,1)
           call sub3(vort(1,2),work1,work2,ntot)
           call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
           call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
           call sub3(vort(1,3),work1,work2,ntot)
      else
           call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
           call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
           call sub3(vort(1,3),work1,work2,ntot)
      endif
      if (ifavg) then
         ifielt = ifield
         ifield = 1
         if (if3d) then
            do idim=1,ldim
               call col2  (vort(1,idim),bm1,ntot)
               call dssum (vort(1,idim),lx1,ly1,lz1)
               call col2  (vort(1,idim),binvm1,ntot)
            enddo
         else
            call col2  (vort(1,3),bm1,ntot)
            call dssum (vort(1,3),lx1,ly1,lz1)
            call col2  (vort(1,3),binvm1,ntot)
         endif
         ifield = ifielt
      endif
      return
      end
      subroutine lorentz_force2(lf,b1,b2,b3)
      include 'SIZE'
      real lf(lx1*ly1*lz1*ldim,lelt)
      real b1(lx1*ly1*lz1,lelt)
      real b2(lx1*ly1*lz1,lelt)
      real b3(lx1*ly1*lz1,lelt)
      integer e
      do e = 1,nelt
         call lorentz_force_e(lf(1,e),b1(1,e),b2(1,e),b3(1,e),e)
      enddo
      return
      end
      subroutine lorentz_force_e(lf,b1,b2,b3,e)
      include 'SIZE'
      include 'DEALIAS'
      include 'GEOM'
      include 'WZ'
      real lf(lx1*ly1*lz1,3)
      real b1(lx1*ly1*lz1)
      real b2(lx1*ly1*lz1)
      real b3(lx1*ly1*lz1)
      integer d,e
      integer icalld
      save    icalld
      data    icalld /0/
      common /ctmp1x/ lfd(lxd*lyd*lzd,3)
     $             ,  bd (lxd*lyd*lzd,3)
     $             ,  cb (lx1*ly1*lz1,3)
     $             ,  cbd(lxd*lyd*lzd,3)
      real lfd

      if (icalld .eq. 0) then
          write(6,*) 'CALL SET PROJ',lx1,lxd
          call setmap (lx1,lxd)
          call setproj(lx1,lxd)
          icalld = icalld + 1
      endif
      call spec_curl_e (cb,b1,b2,b3
     $                 , rxm1(1,1,1,e),rym1(1,1,1,e),rzm1(1,1,1,e)
     $                 , sxm1(1,1,1,e),sym1(1,1,1,e),szm1(1,1,1,e)
     $                 , txm1(1,1,1,e),tym1(1,1,1,e),tzm1(1,1,1,e) )
      do d=1,ldim
         call specmp(cbd(1,d),lxd,cb(1,d),lx1,im1d,im1dt,bd)
      enddo
      call specmp(bd(1,1),lxd,b1,lx1,im1d,im1dt,lfd)
      call specmp(bd(1,2),lxd,b2,lx1,im1d,im1dt,lfd)
      call specmp(bd(1,3),lxd,b3,lx1,im1d,im1dt,lfd)
      nxyzd = lxd*lyd*lzd
      do i=1,nxyzd
         lfd(i,1) = cbd(i,2)*bd(i,3) - cbd(i,3)*bd(i,2)
         lfd(i,2) = cbd(i,3)*bd(i,1) - cbd(i,1)*bd(i,3)
         lfd(i,3) = cbd(i,1)*bd(i,2) - cbd(i,2)*bd(i,1)
      enddo
      call specmp(lf(1,1),lx1,lfd(1,1),lxd,pmd1,pmd1t,cbd)
      call specmp(lf(1,2),lx1,lfd(1,2),lxd,pmd1,pmd1t,cbd)
      call specmp(lf(1,3),lx1,lfd(1,3),lxd,pmd1,pmd1t,cbd)
      nxyz = lx1*ly1*lz1
      do i=1,nxyz
         scale = 1./jacm1(i,1,1,e)
         lf(i,1) = scale*lf(i,1)
         lf(i,2) = scale*lf(i,2)
         lf(i,3) = scale*lf(i,3)
      enddo
      return
      end
      subroutine spec_curl_e (cb,b1,b2,b3,rx,ry,rz,sx,sy,sz,tx,ty,tz)
      include 'SIZE'
      include 'DXYZ'
      real cb(lx1*ly1*lz1,3)
      real b1(1),b2(1),b3(1)
      real rx(1),ry(1),rz(1)
      real sx(1),sy(1),sz(1)
      real tx(1),ty(1),tz(1)
      common /ctmp0x/ br(lx1*ly1*lz1),bs(lx1*ly1*lz1),bt(lx1*ly1*lz1)
      nxyz = lx1*ly1*lz1
      N=lx1-1
      call local_grad3(br,bs,bt,b1,N,1,dxm1,dxtm1)
      do i=1,nxyz
         cb(i,2) =  (br(i)*rz(i)+bs(i)*sz(i)+bt(i)*tz(i))
         cb(i,3) = -(br(i)*ry(i)+bs(i)*sy(i)+bt(i)*ty(i))
      enddo
      call local_grad3(br,bs,bt,b2,N,1,dxm1,dxtm1)
      do i=1,nxyz
         cb(i,1) = -(br(i)*rz(i)+bs(i)*sz(i)+bt(i)*tz(i))
         cb(i,3) =  (br(i)*rx(i)+bs(i)*sx(i)+bt(i)*tx(i))
     $           +  cb(i,3)
      enddo
      call local_grad3(br,bs,bt,b3,N,1,dxm1,dxtm1)
      do i=1,nxyz
         cb(i,1) =  (br(i)*ry(i)+bs(i)*sy(i)+bt(i)*ty(i))
     $           +  cb(i,1)
         cb(i,2) = -(br(i)*rx(i)+bs(i)*sx(i)+bt(i)*tx(i))
     $           +  cb(i,2)
      enddo
      return
      end
      subroutine specx(b,nb,a,na,ba,ab,w)
      include 'SIZE'
      include 'INPUT'
      real b(1),a(1)
      real w(1)
      n=na*na*na
      do i=1,n
         b(i) = a(i)
      enddo
      return
      end
      subroutine phys_to_elsasser(u1,u2,u3,b1,b2,b3,n)
      real u1(1),u2(1),u3(1),b1(1),b2(1),b3(1)
      do i=1,n
         zpx = u1(i) + b1(i)
         zpy = u2(i) + b2(i)
         zpz = u3(i) + b3(i)
         zmx = u1(i) - b1(i)
         zmy = u2(i) - b2(i)
         zmz = u3(i) - b3(i)
         u1(i) = zpx
         u2(i) = zpy
         u3(i) = zpz
         b1(i) = zmx
         b2(i) = zmy
         b3(i) = zmz
      enddo
      return
      end
      subroutine elsasser_to_phys(u1,u2,u3,b1,b2,b3,n)
      real u1(1),u2(1),u3(1),b1(1),b2(1),b3(1)
      do i=1,n
         zpx = 0.5*( u1(i) + b1(i) )
         zpy = 0.5*( u2(i) + b2(i) )
         zpz = 0.5*( u3(i) + b3(i) )
         zmx = 0.5*( u1(i) - b1(i) )
         zmy = 0.5*( u2(i) - b2(i) )
         zmz = 0.5*( u3(i) - b3(i) )
         u1(i) = zpx
         u2(i) = zpy
         u3(i) = zpz
         b1(i) = zmx
         b2(i) = zmy
         b3(i) = zmz
      enddo
      return
      end
      subroutine phys_to_elsasser2(u1,b1,n)
      real u1(1),b1(1)
      do i=1,n
         zpx = u1(i) + b1(i)
         zmx = u1(i) - b1(i)
         u1(i) = zpx
         b1(i) = zmx
      enddo
      return
      end
      subroutine elsasser_to_phys2(u1,b1,n)
      real u1(1),b1(1)
      do i=1,n
         zpx = 0.5*( u1(i) + b1(i) )
         zmx = 0.5*( u1(i) - b1(i) )
         u1(i) = zpx
         b1(i) = zmx
      enddo
      return
      end
      subroutine elsasserh(igeom)
      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'
      include 'GEOM'
      common /scrnt/  besv1 (lbx1,lby1,lbz1,lbelv)
     $ ,              besv2 (lbx1,lby1,lbz1,lbelv)
     $ ,              besv3 (lbx1,lby1,lbz1,lbelv)
      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)
      n  = lx1*ly1*lz1*nelv
      intype = -1
      ifield = 1
      call sethlm   (h1,h2,intype)
      call cresvibp (resv1,resv2,resv3,h1,h2)
      ifield = ifldmhd
      call sethlm   (h1,h2,intype)
      call cresvib  (besv1,besv2,besv3,h1,h2)
      ifield = 1
      call sethlm   (h1,h2,intype)

      call ophinv   (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)

      call opadd2   (vx,vy,vz,dv1,dv2,dv3)

      if (filterType.eq.1) then
         alpha_filt=param(103)
         call q_filter(alpha_filt)
      endif

      call incomprn (vx,vy,vz,pr)

      ifield = ifldmhd
      call sethlm   (h1,h2,intype)

      call ophinv   (dv1,dv2,dv3,besv1,besv2,besv3,h1,h2,tolhv,nmxv)
      call opadd2   (bx,by,bz,dv1,dv2,dv3)

      call incomprn (bx,by,bz,pm)

      return
      end
      subroutine compute_cfl(cfl,u,v,w,dt)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'WZ'
      include 'SOLN'
      real u(lx1,ly1,lz1,nelv),v(lx1,ly1,lz1,nelv),w(lx1,ly1,lz1,nelv)
      common /cfldx/ dri(lx1),dsi(ly1),dti(lz1)
      integer e
      integer icalld
      save    icalld
      data    icalld /0/
      if (icalld.eq.0) then
         icalld=1
         call getdr(dri,zgm1(1,1),lx1)
         call getdr(dsi,zgm1(1,2),ly1)
         if (if3d) call getdr(dti,zgm1(1,3),lz1)
      endif

      cfl = 0.
      l   = 0

      if (if3d) then
         nxyz = lx1*ly1*lz1
         do e=1,nelv
            do k=1,lz1
            do j=1,ly1
            do i=1,lx1
               l = l+1
               ur = ( u(i,j,k,e)*rxm1(i,j,k,e)
     $            +   v(i,j,k,e)*rym1(i,j,k,e)
     $            +   w(i,j,k,e)*rzm1(i,j,k,e) ) * jacmi(l,1)
               us = ( u(i,j,k,e)*sxm1(i,j,k,e)
     $            +   v(i,j,k,e)*sym1(i,j,k,e)
     $            +   w(i,j,k,e)*szm1(i,j,k,e) ) * jacmi(l,1)
               ut = ( u(i,j,k,e)*txm1(i,j,k,e)
     $            +   v(i,j,k,e)*tym1(i,j,k,e)
     $            +   w(i,j,k,e)*tzm1(i,j,k,e) ) * jacmi(l,1)

               cflr = abs(dt*ur*dri(i))
               cfls = abs(dt*us*dsi(j))
               cflt = abs(dt*ut*dti(k))

               cflm = cflr + cfls + cflt
               cfl  = max(cfl,cflm)

               cflf(i,j,k,e) = cflm

            enddo
            enddo
            enddo
         enddo
      else
         nxyz = lx1*ly1
         do e=1,nelv
            do j=1,ly1
            do i=1,lx1
               l = l+1
               ur = ( u(i,j,1,e)*rxm1(i,j,1,e)
     $            +   v(i,j,1,e)*rym1(i,j,1,e) ) * jacmi(l,1)
               us = ( u(i,j,1,e)*sxm1(i,j,1,e)
     $            +   v(i,j,1,e)*sym1(i,j,1,e) ) * jacmi(l,1)

               cflr = abs(dt*ur*dri(i))
               cfls = abs(dt*us*dsi(j))

               cflm = cflr + cfls
               cfl  = max(cfl,cflm)

               cflf(i,j,1,e) = cflm

            enddo
            enddo
         enddo
      endif
      cfl = glmax(cfl,1)
      return
      end
      subroutine getdr(dri,zgm1,lx1)
      real dri(lx1),zgm1(lx1)
      dri(1) = zgm1(2) - zgm1(1)
      do i=2,lx1-1
         dri(i) = 0.5*( zgm1(i+1) - zgm1(i-1) )
      enddo
      dri(lx1) = zgm1(lx1) - zgm1(lx1-1)
      call invcol1(dri,lx1)
      return
      end
      subroutine ophinv(o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'VPROJ'
      include 'TSTEP'

      real o1 (lx1,ly1,lz1,1) , o2 (lx1,ly1,lz1,1) , o3 (lx1,ly1,lz1,1)
      real i1 (lx1,ly1,lz1,1) , i2 (lx1,ly1,lz1,1) , i3 (lx1,ly1,lz1,1)
      real h1 (lx1,ly1,lz1,1) , h2 (lx1,ly1,lz1,1)

      mtmp = param(93)
      do i=1,2*ldim
         ivproj(1,i) = min(mxprev,mtmp) - 1
      enddo

      imesh = 1

      if (ifstrs) then
         matmod = 0
         call hmhzsf  ('NOMG',o1,o2,o3,i1,i2,i3,h1,h2,
     $                  v1mask,v2mask,v3mask,vmult,
     $                  tolh,nmxhi,matmod)
      else
         if (ifield.eq.1) then
            call hsolve ('VELX',o1,i1,h1,h2,v1mask,vmult
     $                         ,imesh,tolh,nmxhi,1
     $                         ,vproj(1,1),ivproj(1,1),binvm1)
            call hsolve ('VELY',o2,i2,h1,h2,v2mask,vmult
     $                         ,imesh,tolh,nmxhi,2
     $                         ,vproj(1,2),ivproj(1,2),binvm1)
            if (if3d)
     $      call hsolve ('VELZ',o3,i3,h1,h2,v3mask,vmult
     $                         ,imesh,tolh,nmxhi,3
     $                         ,vproj(1,3),ivproj(1,3),binvm1)
         elseif (ifield.eq.ifldmhd) then
            call hsolve (' BX ',o1,i1,h1,h2,b1mask,vmult
     $                         ,imesh,tolh,nmxhi,1
     $                         ,vproj(1,4),ivproj(1,4),binvm1)
            call hsolve (' BY ',o2,i2,h1,h2,b2mask,vmult
     $                         ,imesh,tolh,nmxhi,2
     $                         ,vproj(1,5),ivproj(1,5),binvm1)
            if (if3d)
     $      call hsolve (' BZ ',o3,i3,h1,h2,b3mask,vmult
     $                         ,imesh,tolh,nmxhi,3
     $                         ,vproj(1,6),ivproj(1,6),binvm1)
         endif
      endif
      return
      end
      subroutine set_ifbcor
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      common  /nekcb/ cb
      character cb*3

      ifbcor = .true.
      ifield = ifldmhd

      nface  = 2*ldim
      do iel=1,nelv
      do ifc=1,nface
         cb = cbc(ifc,iel,ifield)
         if  (cb.eq.'ndd' .or. cb.eq.'dnd' .or. cb.eq.'ddn')
     $           ifbcor = .false.
      enddo
      enddo

      call gllog(ifbcor , .false.)

      if (nio.eq.0)  write (6,*) 'IFBCOR   =',ifbcor

      return
      end
      subroutine setrhsp(p,h1,h2,h2inv,pset,niprev)
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
      real pset (lx2*ly2*lz2*lelv,mxprev)

      parameter (ltot2=lx2*ly2*lz2*lelv)
      common /orthox/ pbar(ltot2),pnew(ltot2)
      common /orthos/ alpha(mxprev),work(mxprev)
      common /orthoi/ nprev,mprev

      nprev = niprev
      if (nprev.eq.0) return

      ntot2  = lx2*ly2*lz2*nelv
      alpha1 = glsc3(p,p,bm2inv,ntot2)
      if (alpha1.gt.0) then
         alpha1 = sqrt(alpha1/volvm2)
      else
         return
      endif

      call updrhse(p,h1,h2,h2inv,ierr)

      do i=1,nprev
         alpha(i) = vlsc2(p,pset(1,i),ntot2)
      enddo
      call gop(alpha,work,'+  ',nprev)

      call rzero(pbar,ntot2)
      do i=1,nprev
         call add2s2(pbar,pset(1,i),alpha(i),ntot2)
      enddo
      intetype = 1
      call cdabdtp(pnew,pbar,h1,h2,h2inv,intetype)
      call sub2   (p,pnew,ntot2)

      alpha2 = glsc3(p,p,bm2inv,ntot2)
      if (alpha2.gt.0) then
         alpha2 = sqrt(alpha2/volvm2)
         ratio  = alpha1/alpha2
         n10=min(10,nprev)
   11    format(2i5,' alpha:',1p10e12.4)
   12    format(i6,i4,1p3e12.4,' alph12')

         if (nio.eq.0) write(6,13) istep,'  Project PRES  ',
     &                         alpha2,alpha1,ratio,nprev,mxprev
   13    format(i11,a,6x,1p3e13.4,i4,i4)
      endif

      return
      end
      subroutine gensolnp(p,h1,h2,h2inv,pset,nprev)
      include 'SIZE'
      include 'INPUT'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
      real pset (lx2*ly2*lz2*lelv,mxprev)

      parameter (ltot2=lx2*ly2*lz2*lelv)
      common /orthox/ pbar(ltot2),pnew(ltot2)
      common /orthos/ alpha(mxprev),work(mxprev)

      mprev=param(93)
      mprev=min(mprev,mxprev)

      ntot2=lx2*ly2*lz2*nelv

      if (nprev.lt.mprev) then
         nprev = nprev+1
         call copy  (pset(1,nprev),p,ntot2)
         call add2  (p,pbar,ntot2)
         call econjp(pset,nprev,h1,h2,h2inv,ierr)

         if (ierr.eq.1) then
          nprev = 1
          call copy  (pset(1,nprev),p,ntot2)
          call econjp(pset,nprev,h1,h2,h2inv,ierr)
        endif

      else
         nprev = 1
         call add2  (p,pbar,ntot2)
         call copy  (pset(1,nprev),p,ntot2)
         call econjp(pset,nprev,h1,h2,h2inv,ierr)
      endif

      return
      end
      subroutine econjp(pset,nprev,h1,h2,h2inv,ierr)

      include 'SIZE'
      include 'INPUT'

      real p    (lx2,ly2,lz2,lelv)
      real h1   (lx1,ly1,lz1,lelv)
      real h2   (lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)
      real pset (lx2*ly2*lz2*lelv,mxprev)

      parameter (ltot2=lx2*ly2*lz2*lelv)
      common /orthox/ pbar(ltot2),pnew(ltot2)
      common /orthos/ alpha(mxprev),work(mxprev)

      ierr  = 0

      ntot2=lx2*ly2*lz2*nelv

      npass=1
      do ipass=1,npass

         intetype=1
         call cdabdtp(pnew,pset(1,nprev),h1,h2,h2inv,intetype)
         alphad = glsc2(pnew,pset(1,nprev),ntot2)

         nprev1 = nprev-1
         do i=1,nprev1
            alpha(i) = vlsc2(pnew,pset(1,i),ntot2)
         enddo
         if (nprev1.gt.0) call gop(alpha,work,'+  ',nprev1)

         do i=1,nprev1
            alpham = -alpha(i)
            call add2s2(pset(1,nprev),pset(1,i),alpham,ntot2)
            alphad = alphad - alpha(i)**2
         enddo
      enddo
      if (alphad.le.0) then
         write(6,*) 'ERROR:  alphad .le. 0 in econjp',alphad,nprev
         ierr = 1
         return
      endif
      alphad = 1./sqrt(alphad)
      call cmult(pset(1,nprev),alphad,ntot2)

      return
      end
      subroutine advab_elsasser_fast
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrns/ wk(2*ltd)
     $             , fx(lxy),fy(lxy),fz(lxy)
     $             , gx(lxy),gy(lxy),gz(lxy)
     $             , zr(ltd),zs(ltd),zt(ltd)
     $             , tr(ltd,3),zp(ltd,3),zm(ltd,3)

      integer e
      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.0) call set_dealias_rx
      icalld = icalld + 1

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      if (if3d) then
         do e=1,nelv

            call add3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,1),wk,lx1,lxd,if3d,0)

            call add3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,2),wk,lx1,lxd,if3d,0)

            call add3(wk,vz(1,1,1,e),bz(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,3),wk,lx1,lxd,if3d,0)

            call sub3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,1),wk,lx1,lxd,if3d,0)

            call sub3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,2),wk,lx1,lxd,if3d,0)

            call sub3(wk,vz(1,1,1,e),bz(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,3),wk,lx1,lxd,if3d,0)

            do i=1,nxyzd
               tr(i,1)=
     $            rx(i,1,e)*zm(i,1)+rx(i,2,e)*zm(i,2)+rx(i,3,e)*zm(i,3)
               tr(i,2)=
     $            rx(i,4,e)*zm(i,1)+rx(i,5,e)*zm(i,2)+rx(i,6,e)*zm(i,3)
               tr(i,3)=
     $            rx(i,7,e)*zm(i,1)+rx(i,8,e)*zm(i,2)+rx(i,9,e)*zm(i,3)
            enddo

            call grad_rst(zr,zs,zt,zp(1,1),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fx,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zp(1,2),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fy,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zp(1,3),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(fz,wk,lx1,lxd,if3d,1)

            do i=1,nxyzd
               tr(i,1)=
     $            rx(i,1,e)*zp(i,1)+rx(i,2,e)*zp(i,2)+rx(i,3,e)*zp(i,3)
               tr(i,2)=
     $            rx(i,4,e)*zp(i,1)+rx(i,5,e)*zp(i,2)+rx(i,6,e)*zp(i,3)
               tr(i,3)=
     $            rx(i,7,e)*zp(i,1)+rx(i,8,e)*zp(i,2)+rx(i,9,e)*zp(i,3)
            enddo

            call grad_rst(zr,zs,zt,zm(1,1),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gx,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zm(1,2),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gy,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zm(1,3),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)+tr(i,3)*zt(i)
            enddo
            call intp_rstd(gz,wk,lx1,lxd,if3d,1)

            tmp = -0.5
            do i=1,nxyz1
               bfx(i,1,1,e) = bfx(i,1,1,e) + tmp*( fx(i) + gx(i) )
               bmx(i,1,1,e) = bmx(i,1,1,e) + tmp*( fx(i) - gx(i) )

               bfy(i,1,1,e) = bfy(i,1,1,e) + tmp*( fy(i) + gy(i) )
               bmy(i,1,1,e) = bmy(i,1,1,e) + tmp*( fy(i) - gy(i) )

               bfz(i,1,1,e) = bfz(i,1,1,e) + tmp*( fz(i) + gz(i) )
               bmz(i,1,1,e) = bmz(i,1,1,e) + tmp*( fz(i) - gz(i) )
            enddo

         enddo

      else
         do e=1,nelv

            call add3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,1),wk,lx1,lxd,if3d,0)

            call add3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zp(1,2),wk,lx1,lxd,if3d,0)

            call sub3(wk,vx(1,1,1,e),bx(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,1),wk,lx1,lxd,if3d,0)

            call sub3(wk,vy(1,1,1,e),by(1,1,1,e),nxyz1)
            call intp_rstd(zm(1,2),wk,lx1,lxd,if3d,0)

            do i=1,nxyzd
               tr(i,1)=
     $            rx(i,1,e)*zm(i,1)+rx(i,2,e)*zm(i,2)
               tr(i,2)=
     $            rx(i,3,e)*zm(i,1)+rx(i,4,e)*zm(i,2)
            enddo

            call grad_rst(zr,zs,zt,zp(1,1),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(fx,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zp(1,2),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(fy,wk,lx1,lxd,if3d,1)

            do i=1,nxyzd
               tr(i,1)=
     $            rx(i,1,e)*zp(i,1)+rx(i,2,e)*zp(i,2)
               tr(i,2)=
     $            rx(i,3,e)*zp(i,1)+rx(i,4,e)*zp(i,2)
            enddo

            call grad_rst(zr,zs,zt,zm(1,1),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(gx,wk,lx1,lxd,if3d,1)

            call grad_rst(zr,zs,zt,zm(1,2),lxd,if3d)
            do i=1,nxyzd
               wk(i) = tr(i,1)*zr(i)+tr(i,2)*zs(i)
            enddo
            call intp_rstd(gy,wk,lx1,lxd,if3d,1)

            tmp = -0.5
            do i=1,nxyz1
               bfx(i,1,1,e) = bfx(i,1,1,e) + tmp*( fx(i) + gx(i) )
               bmx(i,1,1,e) = bmx(i,1,1,e) + tmp*( fx(i) - gx(i) )

               bfy(i,1,1,e) = bfy(i,1,1,e) + tmp*( fy(i) + gy(i) )
               bmy(i,1,1,e) = bmy(i,1,1,e) + tmp*( fy(i) - gy(i) )
            enddo
         enddo
      endif

      return
      end
      subroutine cfl_check
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
     $ ,             tb1 (lx1,ly1,lz1,lelv)
     $ ,             tb2 (lx1,ly1,lz1,lelv)
     $ ,             tb3 (lx1,ly1,lz1,lelv)

      call opcopy(ta1,ta2,ta3,vx,vy,vz)
      call opcopy(tb1,tb2,tb3,bx,by,bz)

      ntot1 = lx1*ly1*lz1*nelv
      call phys_to_elsasser(ta1,ta2,ta3,tb1,tb2,tb3,ntot1)

      call compute_cfl(cflp,ta1,ta2,ta3,dt)
      call compute_cfl(cflm,tb1,tb2,tb3,dt)

      courno = max(cflp,cflm)

      if (nio.eq.0) write(6,1) istep,time,dt,cflp,cflm
    1 format(i9,1p4e15.7,' CFL')

      return
      end
