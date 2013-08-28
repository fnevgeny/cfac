      program dummy
      implicit none

c     DB filename
      character (len=256) dbname

c     species
      integer anum, nele_min, nele_max
      double precision mass

      integer ndim, rtdim, aidim, cedim, cidim, pidim

      integer ierr

c     Sink subroutines for handling (storing) data      
      external l_sink, rt_sink, ai_sink, ct_sink

      dbname = 'Ge.db'
      
      nele_min = 1
      nele_max = 2

c     Initialization; obtain various dimensions for dynamic allocation
      call cfacdb_init(dbname, nele_min, nele_max,
     &                 ndim, rtdim, aidim, cedim, cidim, pidim,
     &                 ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_init() failed with ierr = ', ierr
          stop
      endif
      
c     Get species properties
      call cfacdb_species(anum, mass, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_species() failed with ierr = ', ierr
          stop
      endif
      
      write(*, *) '----------------------------------------------'
      
      write(*, 901) anum, mass
      write(*, 902) ndim, rtdim, aidim, cedim, cidim, pidim
      
      write(*, *) '----------------------------------------------'
      write(*, *) '   #   id    E nele    g  vn  vl   p      name'
      write(*, *) '----------------------------------------------'

      call cfacdb_levels(l_sink, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_levels() failed with ierr = ', ierr
          stop
      endif
      
      write(*, *) '----------------------------------------------'

c     Get radiative transitions
      call cfacdb_rtrans(rt_sink, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_rtrans() failed with ierr = ', ierr
          stop
      endif
      
c     Get AI transitions
      call cfacdb_aitrans(ai_sink, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_aitrans() failed with ierr = ', ierr
          stop
      endif
      
c     Get collisional (CE, CI, PI) transitions
      call cfacdb_ctrans(ct_sink, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_ctrans() failed with ierr = ', ierr
          stop
      endif
      
c     Close the DB and free associated structures
      call cfacdb_close()
      
c---------
      
      stop

 901  format(' Species: anum =', i3, '; mass =', f7.2)
 902  format(' ndim =', i5, '; rtdim =', i5, '; aidim =', i5,
     &       '; cedim =', i5, '; cidim =', i5,'; pidim =', i5)

      end


      subroutine l_sink(id, e, nele, g, vn, vl, p, name, ncmplx, sname)
      implicit none
      integer id, nele, g, vn, vl, p
      double precision e
      character (len=*) name, ncmplx, sname

      write (*,920) id, e, nele, g, vn, vl, p, name
 920  format(i5, f10.3, i5, i5, i4, i4, i4, ' ', a20)

      return
      end

      subroutine rt_sink(i, j, mpole, gf)
      implicit none
      integer i, j, mpole
      double precision gf

      write (*,921) i, j, mpole, gf
 921  format(i5, ' -> ', i5, ' mpole = ', i2, ', gf = ', g10.3)

      return
      end

      subroutine ai_sink(i, j, rate)
      implicit none
      integer i, j
      double precision rate

      write (*,922) i, j, rate
 922  format(i5, ' -> ', i5, ', ai = ', g10.3)

      return
      end

      subroutine ct_sink(i, j, ctype, ap0, ap1, nd, e, d)
      implicit none
      integer i, j, ctype, nd
      double precision e(*), d(*), ap0, ap1

      write (*,923) i, j, ctype, ap0, ap1, e(1:nd), d(1:nd)
 923  format(i5, ' -> ', i5, ', type = ', i1, ' ap0,1 = [',
     &       g10.3, g10.3, ']', 12g10.3)

      return
      end
