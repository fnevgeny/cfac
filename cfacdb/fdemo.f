c A F77 demo showing use of the CFACDB API
c 
c Copyright (C) 2013 Evgeny Stambulchik
c 
c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
c
      program fdemo
      implicit none

c     DB filename
      character (len=256) dbname

c     species
      integer anum, nele_min, nele_max
      double precision mass

      integer ndim, rtdim, aidim, cedim, cidim, pidim
      
      double precision T

      integer ierr

c     Sink subroutines for handling (storing) data      
      external l_sink, rt_sink, ai_sink, ct_sink, cr_sink

      dbname = 'Ge.db'
      
      nele_min = 0
      nele_max = 3

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
      write(*, *) '   #         E nele    g  vn  vl   p      name'
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
c      call cfacdb_ctrans(ct_sink, ierr)

c     Get collisional (CE, CI, PI) rates
c    (NB: the statweight of the ini state is included!!!)
      T = 2000.0/27.211
      call cfacdb_crates(T, cr_sink, ierr)
      if (ierr .ne. 0) then
          print *, 'cfacdb_crates() failed with ierr = ', ierr
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

      subroutine rt_sink(i, j, mpole, gf, uta_de, uta_sd)
      implicit none
      integer i, j, mpole
      double precision gf, uta_de, uta_sd

      write (*,921) i, j, mpole, gf, uta_de, uta_sd
 921  format(i5, ' -> ', i5, ' mpole = ', i2, ', gf = ', g10.3,
     c       ', uta_de = ', g10.3, ', uta_sd = ', g10.3)

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

      subroutine ct_sink(i, j, ctype, kl, ap0, ap1, ap2, ap3, nd, e, d)
      implicit none
      integer i, j, ctype, kl, nd
      double precision e(*), d(*), ap0, ap1, ap2, ap3

      write (*,923) i, j, ctype, kl, ap0, ap1, ap2, ap3,
     &              e(1:nd), d(1:nd)
 923  format(i5, ' -> ', i5, ', type = ', i1, ', kl = ', i1,
     &       ' ap0-3 = [', 4g10.3, ']', 12g10.3)

      return
      end

      subroutine cr_sink(i, j, ctype, ratec)
      implicit none
      integer i, j, ctype
      double precision ratec

      write (*,924) i, j, ctype, ratec
 924  format(i5, ' -> ', i5, ', type = ', i1, ' ratec = ', g10.3)

      return
      end
