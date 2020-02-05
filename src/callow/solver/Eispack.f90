!------------------------------------------------------------------------------!
!> @file  Eispack.f90
!> @brief Fortran 90 conversion of several Eispack routines
!>
!> Routines are included to compute eigenvalues and vectors of standard and
!> generalized eigenvalue problems.  All subroutines have their original
!> comments.
!>
!> Note, the subroutines RG and RGG drive the standard and generalized
!> problem solves, respectively.
!------------------------------------------------------------------------------!

      SUBROUTINE balanc (nm, n, a, low, igh, scale)
!
      INTEGER i, j, k, l, m, n, jj, nm, igh, low, iexc
      DOUBLE PRECISION a (nm, n), scale (n)
      DOUBLE PRECISION c, f, g, r, s, b2, radix
      LOGICAL noconv
!
!     this subroutine is a translation of the algol procedure balance,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine balances a real matrix and isolates
!     eigenvalues whenever possible.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the input matrix to be balanced.
!
!     on output
!
!        a contains the balanced matrix.
!
!        low and igh are two integers such that a(i,j)
!          is equal to zero if
!           (1) i is greater than j and
!           (2) j=1,...,low-1 or i=igh+1,...,n.
!
!        scale contains information determining the
!           permutations and scaling factors used.
!
!     suppose that the principal submatrix in rows low through igh
!     has been balanced, that p(j) denotes the index interchanged
!     with j during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by d(i,j).  then
!        scale(j) = p(j),    for j = 1,...,low-1
!                 = d(j,j),      j = low,...,igh
!                 = p(j)         j = igh+1,...,n.
!     the order in which the interchanges are made is n to igh+1,
!     then 1 to low-1.
!
!     note that 1 is returned for igh if igh is zero formally.
!
!     the algol procedure exc contained in balance appears in
!     balanc  in line.  (note that the algol roles of identifiers
!     k,l have been reversed.)
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      radix = 16.0d0
!
      b2 = radix * radix
      k = 1
      l = n
      GOTO 100
!     .......... in-line procedure for row and
!                column exchange ..........
   20 scale (m) = j
      IF (j.eq.m) goto 50
!
      DO 30 i = 1, l
        f = a (i, j)
        a (i, j) = a (i, m)
        a (i, m) = f
   30 END DO
!
      DO 40 i = k, n
        f = a (j, i)
        a (j, i) = a (m, i)
        a (m, i) = f
   40 END DO
!
   50 GOTO (80, 130), iexc
!     .......... search for rows isolating an eigenvalue
!                and push them down ..........
   80 IF (l.eq.1) goto 280
      l = l - 1
!     .......... for j=l step -1 until 1 do -- ..........
  100 DO 120 jj = 1, l
        j = l + 1 - jj
!
        DO 110 i = 1, l
          IF (i.eq.j) goto 110
          IF (a (j, i) .ne.0.0d0) goto 120
  110   END DO
!
        m = l
        iexc = 1
        GOTO 20
  120 END DO
!
      GOTO 140
!     .......... search for columns isolating an eigenvalue
!                and push them left ..........
  130 k = k + 1
!
  140 DO 170 j = k, l
!
        DO 150 i = k, l
          IF (i.eq.j) goto 150
          IF (a (i, j) .ne.0.0d0) goto 170
  150   END DO
!
        m = k
        iexc = 2
        GOTO 20
  170 END DO
!     .......... now balance the submatrix in rows k to l ..........
      DO 180 i = k, l
  180 scale (i) = 1.0d0
!     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
!
      DO 270 i = k, l
        c = 0.0d0
        r = 0.0d0
!
        DO 200 j = k, l
          IF (j.eq.i) goto 200
          c = c + dabs (a (j, i) )
          r = r + dabs (a (i, j) )
  200   END DO
!     .......... guard against zero c or r due to underflow ..........
        IF (c.eq.0.0d0.or.r.eq.0.0d0) goto 270
        g = r / radix
        f = 1.0d0
        s = c + r
  210   IF (c.ge.g) goto 220
        f = f * radix
        c = c * b2
        GOTO 210
  220   g = r * radix
  230   IF (c.lt.g) goto 240
        f = f / radix
        c = c / b2
        GOTO 230
!     .......... now balance ..........
  240   IF ( (c + r) / f.ge.0.95d0 * s) goto 270
        g = 1.0d0 / f
        scale (i) = scale (i) * f
        noconv = .true.
!
        DO 250 j = k, n
  250   a (i, j) = a (i, j) * g
!
        DO 260 j = 1, l
  260   a (j, i) = a (j, i) * f
!
  270 END DO
!
      IF (noconv) goto 190
!
  280 low = k
      igh = l
      RETURN
      END SUBROUTINE balanc

      SUBROUTINE balbak (nm, n, low, igh, scale, m, z)
!
      INTEGER i, j, k, m, n, ii, nm, igh, low
      DOUBLE PRECISION scale (n), z (nm, m)
      DOUBLE PRECISION s
!
!     this subroutine is a translation of the algol procedure balbak,
!     num. math. 13, 293-304(1969) by parlett and reinsch.
!     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
!
!     this subroutine forms the eigenvectors of a real general
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  balanc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by  balanc.
!
!        scale contains information determining the permutations
!          and scaling factors used by  balanc.
!
!        m is the number of columns of z to be back transformed.
!
!        z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first m columns.
!
!     on output
!
!        z contains the real and imaginary parts of the
!          transformed eigenvectors in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      IF (m.eq.0) goto 200
      IF (igh.eq.low) goto 120
!
      DO 110 i = low, igh
        s = scale (i)
!     .......... left hand eigenvectors are back transformed
!                if the foregoing statement is replaced by
!                s=1.0d0/scale(i). ..........
        DO 100 j = 1, m
  100   z (i, j) = z (i, j) * s
!
  110 END DO
!     ......... for i=low-1 step -1 until 1,
!               igh+1 step 1 until n do -- ..........
  120 DO 140 ii = 1, n
        i = ii
        IF (i.ge.low.and.i.le.igh) goto 140
        IF (i.lt.low) i = low - ii
        k = scale (i)
        IF (k.eq.i) goto 140
!
        DO 130 j = 1, m
          s = z (i, j)
          z (i, j) = z (k, j)
          z (k, j) = s
  130   END DO
!
  140 END DO
!
  200 RETURN
      END SUBROUTINE balbak

      SUBROUTINE elmhes (nm, n, low, igh, a, int)
!
      INTEGER i, j, m, n, la, nm, igh, kp1, low, mm1, mp1
      DOUBLE PRECISION a (nm, n)
      DOUBLE PRECISION x, y
      INTEGER int (igh)
!
!     this subroutine is a translation of the algol procedure elmhes,
!     num. math. 12, 349-368(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!     given a real general matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     low through igh to upper hessenberg form by
!     stabilized elementary similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the input matrix.
!
!     on output
!
!        a contains the hessenberg matrix.  the multipliers
!          which were used in the reduction are stored in the
!          remaining triangle under the hessenberg matrix.
!
!        int contains information on the rows and columns
!          interchanged in the reduction.
!          only elements low through igh are used.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      IF (la.lt.kp1) goto 200
!
      DO 180 m = kp1, la
        mm1 = m - 1
        x = 0.0d0
        i = m
!
        DO 100 j = m, igh
          IF (dabs (a (j, mm1) ) .le.dabs (x) ) goto 100
          x = a (j, mm1)
          i = j
  100   END DO
!
        int (m) = i
        IF (i.eq.m) goto 130
!     .......... interchange rows and columns of a ..........
        DO 110 j = mm1, n
          y = a (i, j)
          a (i, j) = a (m, j)
          a (m, j) = y
  110   END DO
!
        DO 120 j = 1, igh
          y = a (j, i)
          a (j, i) = a (j, m)
          a (j, m) = y
  120   END DO
!     .......... end interchange ..........
  130   IF (x.eq.0.0d0) goto 180
        mp1 = m + 1
!
        DO 160 i = mp1, igh
          y = a (i, mm1)
          IF (y.eq.0.0d0) goto 160
          y = y / x
          a (i, mm1) = y
!
          DO 140 j = m, n
  140     a (i, j) = a (i, j) - y * a (m, j)
!
          DO 150 j = 1, igh
  150     a (j, m) = a (j, m) + y * a (j, i)
!
  160   END DO
!
  180 END DO
!
  200 RETURN
      END SUBROUTINE elmhes

      SUBROUTINE eltran (nm, n, low, igh, a, int, z)
!
      INTEGER i, j, n, kl, mm, mp, nm, igh, low, mp1
      DOUBLE PRECISION a (nm, igh), z (nm, n)
      INTEGER int (igh)
!
!     this subroutine is a translation of the algol procedure elmtrans,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     real general matrix to upper hessenberg form by  elmhes.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        a contains the multipliers which were used in the
!          reduction by  elmhes  in its lower triangle
!          below the subdiagonal.
!
!        int contains information on the rows and columns
!          interchanged in the reduction by  elmhes.
!          only elements low through igh are used.
!
!     on output
!
!        z contains the transformation matrix produced in the
!          reduction by  elmhes.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z to identity matrix ..........
      DO 80 j = 1, n
!
        DO 60 i = 1, n
   60   z (i, j) = 0.0d0
!
        z (j, j) = 1.0d0
   80 END DO
!
      kl = igh - low - 1
      IF (kl.lt.1) goto 200
!     .......... for mp=igh-1 step -1 until low+1 do -- ..........
      DO 140 mm = 1, kl
        mp = igh - mm
        mp1 = mp + 1
!
        DO 100 i = mp1, igh
  100   z (i, mp) = a (i, mp - 1)
!
        i = int (mp)
        IF (i.eq.mp) goto 140
!
        DO 130 j = mp, igh
          z (mp, j) = z (i, j)
          z (i, j) = 0.0d0
  130   END DO
!
        z (i, mp) = 1.0d0
  140 END DO
!
  200 RETURN
      END SUBROUTINE eltran

      SUBROUTINE hqr (nm, n, low, igh, h, wr, wi, ierr)
!  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
!
      INTEGER i, j, k, l, m, n, en, ll, mm, na, nm, igh, itn, its, low, &
      mp2, enm2, ierr
      DOUBLE PRECISION h (nm, n), wr (n), wi (n)
      DOUBLE PRECISION p, q, r, s, t, w, x, y, zz, norm, tst1, tst2
      LOGICAL notlas
!
!     this subroutine is a translation of the algol procedure hqr,
!     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!     this subroutine finds the eigenvalues of a real
!     upper hessenberg matrix by the qr method.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.  information about
!          the transformations used in the reduction to hessenberg
!          form by  elmhes  or  orthes, if performed, is stored
!          in the remaining triangle under the hessenberg matrix.
!
!     on output
!
!        h has been destroyed.  therefore, it must be saved
!          before calling  hqr  if subsequent calculation and
!          back transformation of eigenvectors is to be performed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated september 1989.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0d0
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      DO 50 i = 1, n
!
        DO 40 j = k, n
   40   norm = norm + dabs (h (i, j) )
!
        k = i
        IF (i.ge.low.and.i.le.igh) goto 50
        wr (i) = h (i, i)
        wi (i) = 0.0d0
   50 END DO
!
      en = igh
      t = 0.0d0
      itn = 30 * n
!     .......... search for next eigenvalues ..........
   60 IF (en.lt.low) goto 1001
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 DO 80 ll = low, en
        l = en + low - ll
        IF (l.eq.low) goto 100
        s = dabs (h (l - 1, l - 1) ) + dabs (h (l, l) )
        IF (s.eq.0.0d0) s = norm
        tst1 = s
        tst2 = tst1 + dabs (h (l, l - 1) )
        IF (tst2.eq.tst1) goto 100
   80 END DO
!     .......... form shift ..........
  100 x = h (en, en)
      IF (l.eq.en) goto 270
      y = h (na, na)
      w = h (en, na) * h (na, en)
      IF (l.eq.na) goto 280
      IF (itn.eq.0) goto 1000
      IF (its.ne.10.and.its.ne.20) goto 130
!     .......... form exceptional shift ..........
      t = t + x
!
      DO 120 i = low, en
  120 h (i, i) = h (i, i) - x
!
      s = dabs (h (en, na) ) + dabs (h (na, enm2) )
      x = 0.75d0 * s
      y = x
      w = - 0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      DO 140 mm = l, enm2
        m = enm2 + l - mm
        zz = h (m, m)
        r = x - zz
        s = y - zz
        p = (r * s - w) / h (m + 1, m) + h (m, m + 1)
        q = h (m + 1, m + 1) - zz - r - s
        r = h (m + 2, m + 1)
        s = dabs (p) + dabs (q) + dabs (r)
        p = p / s
        q = q / s
        r = r / s
        IF (m.eq.l) goto 150
        tst1 = dabs (p) * (dabs (h (m - 1, m - 1) ) + dabs (zz) + dabs (&
        h (m + 1, m + 1) ) )
        tst2 = tst1 + dabs (h (m, m - 1) ) * (dabs (q) + dabs (r) )
        IF (tst2.eq.tst1) goto 150
  140 END DO
!
  150 mp2 = m + 2
!
      DO 160 i = mp2, en
        h (i, i - 2) = 0.0d0
        IF (i.eq.mp2) goto 160
        h (i, i - 3) = 0.0d0
  160 END DO
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      DO 260 k = m, na
        notlas = k.ne.na
        IF (k.eq.m) goto 170
        p = h (k, k - 1)
        q = h (k + 1, k - 1)
        r = 0.0d0
        IF (notlas) r = h (k + 2, k - 1)
        x = dabs (p) + dabs (q) + dabs (r)
        IF (x.eq.0.0d0) goto 260
        p = p / x
        q = q / x
        r = r / x
  170   s = dsign (dsqrt (p * p + q * q + r * r), p)
        IF (k.eq.m) goto 180
        h (k, k - 1) = - s * x
        GOTO 190
  180   IF (l.ne.m) h (k, k - 1) = - h (k, k - 1)
  190   p = p + s
        x = p / s
        y = q / s
        zz = r / s
        q = q / p
        r = r / p
        IF (notlas) goto 225
!     .......... row modification ..........
        DO 200 j = k, EN
          p = h (k, j) + q * h (k + 1, j)
          h (k, j) = h (k, j) - p * x
          h (k + 1, j) = h (k + 1, j) - p * y
  200   END DO
!
        j = min0 (en, k + 3)
!     .......... column modification ..........
        DO 210 i = L, j
          p = x * h (i, k) + y * h (i, k + 1)
          h (i, k) = h (i, k) - p
          h (i, k + 1) = h (i, k + 1) - p * q
  210   END DO
        GOTO 255
  225   CONTINUE
!     .......... row modification ..........
        DO 230 j = k, EN
          p = h (k, j) + q * h (k + 1, j) + r * h (k + 2, j)
          h (k, j) = h (k, j) - p * x
          h (k + 1, j) = h (k + 1, j) - p * y
          h (k + 2, j) = h (k + 2, j) - p * zz
  230   END DO
!
        j = min0 (en, k + 3)
!     .......... column modification ..........
        DO 240 i = L, j
          p = x * h (i, k) + y * h (i, k + 1) + zz * h (i, k + 2)
          h (i, k) = h (i, k) - p
          h (i, k + 1) = h (i, k + 1) - p * q
          h (i, k + 2) = h (i, k + 2) - p * r
  240   END DO
  255   CONTINUE
!
  260 END DO
!
      GOTO 70
!     .......... one root found ..........
  270 wr (en) = x + t
      wi (en) = 0.0d0
      en = na
      GOTO 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt (dabs (q) )
      x = x + t
      IF (q.lt.0.0d0) goto 320
!     .......... real pair ..........
      zz = p + dsign (zz, p)
      wr (na) = x + zz
      wr (en) = wr (na)
      IF (zz.ne.0.0d0) wr (en) = x - w / zz
      wi (na) = 0.0d0
      wi (en) = 0.0d0
      GOTO 330
!     .......... complex pair ..........
  320 wr (na) = x + p
      wr (en) = x + p
      wi (na) = zz
      wi (en) = - zz
  330 en = enm2
      GOTO 60
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 RETURN
      END SUBROUTINE hqr

      SUBROUTINE rg (nm, n, a, wr, wi, matz, z, iv1, fv1, ierr)
!
      INTEGER n, nm, is1, is2, ierr, matz
      DOUBLE PRECISION a (nm, n), wr (n), wi (n), z (nm, n), fv1 (n)
      INTEGER iv1 (n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real general matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        wr  and  wi  contain the real and imaginary parts,
!        respectively, of the eigenvalues.  complex conjugate
!        pairs of eigenvalues appear consecutively with the
!        eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for hqr
!           and hqr2.  the normal completion code is zero.
!
!        iv1  and  fv1  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      IF (n.le.nm) goto 10
      ierr = 10 * n
      GOTO 50
!
   10 CALL balanc (nm, n, a, is1, is2, fv1)
      CALL elmhes (nm, n, is1, is2, a, iv1)
      IF (matz.ne.0) goto 20
!     .......... find eigenvalues only ..........
      CALL hqr (nm, n, is1, is2, a, wr, wi, ierr)
      GOTO 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 CALL eltran (nm, n, is1, is2, a, iv1, z)
      CALL hqr2 (nm, n, is1, is2, a, wr, wi, z, ierr)
      IF (ierr.ne.0) goto 50
      CALL balbak (nm, n, is1, is2, fv1, n, z)
   50 RETURN
      END SUBROUTINE rg

      SUBROUTINE qzhes (nm, n, a, b, matz, z)
!
      INTEGER i, j, k, l, n, lb, l1, nm, nk1, nm1, nm2
      DOUBLE PRECISION a (nm, n), b (nm, n), z (nm, n)
      DOUBLE PRECISION r, s, t, u1, u2, v1, v2, rho
      LOGICAL matz
!
!     this subroutine is the first step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real general matrices and
!     reduces one of them to upper hessenberg form and the other
!     to upper triangular form using orthogonal transformations.
!     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real general matrix.
!
!        b contains a real general matrix.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!     on output
!
!        a has been reduced to upper hessenberg form.  the elements
!          below the first subdiagonal have been set to zero.
!
!        b has been reduced to upper triangular form.  the elements
!          below the main diagonal have been set to zero.
!
!        z contains the product of the right hand transformations if
!          matz has been set to .true.  otherwise, z is not referenced.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
!     .......... initialize z ..........
      IF (.not.matz) goto 10
!
      DO 3 j = 1, n
!
        DO 2 i = 1, n
          z (i, j) = 0.0d0
    2   END DO
!
        z (j, j) = 1.0d0
    3 END DO
!     .......... reduce b to upper triangular form ..........
   10 IF (n.le.1) goto 170
      nm1 = n - 1
!
      DO 100 l = 1, nm1
        l1 = l + 1
        s = 0.0d0
!
        DO 20 i = l1, n
          s = s + dabs (b (i, l) )
   20   END DO
!
        IF (s.eq.0.0d0) goto 100
        s = s + dabs (b (l, l) )
        r = 0.0d0
!
        DO 25 i = l, n
          b (i, l) = b (i, l) / s
          r = r + b (i, l) **2
   25   END DO
!
        r = dsign (dsqrt (r), b (l, l) )
        b (l, l) = b (l, l) + r
        rho = r * b (l, l)
!
        DO 50 j = l1, n
          t = 0.0d0
!
          DO 30 i = l, n
            t = t + b (i, l) * b (i, j)
   30     END DO
!
          t = - t / rho
!
          DO 40 i = l, n
            b (i, j) = b (i, j) + t * b (i, l)
   40     END DO
!
   50   END DO
!
        DO 80 j = 1, n
          t = 0.0d0
!
          DO 60 i = l, n
            t = t + b (i, l) * a (i, j)
   60     END DO
!
          t = - t / rho
!
          DO 70 i = l, n
            a (i, j) = a (i, j) + t * b (i, l)
   70     END DO
!
   80   END DO
!
        b (l, l) = - s * r
!
        DO 90 i = l1, n
          b (i, l) = 0.0d0
   90   END DO
!
  100 END DO
!     .......... reduce a to upper hessenberg form, while
!                keeping b triangular ..........
      IF (n.eq.2) goto 170
      nm2 = n - 2
!
      DO 160 k = 1, nm2
        nk1 = nm1 - k
!     .......... for l=n-1 step -1 until k+1 do -- ..........
        DO 150 lb = 1, nk1
          l = n - lb
          l1 = l + 1
!     .......... zero a(l+1,k) ..........
          s = dabs (a (l, k) ) + dabs (a (l1, k) )
          IF (s.eq.0.0d0) goto 150
          u1 = a (l, k) / s
          u2 = a (l1, k) / s
          r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
          v1 = - (u1 + r) / r
          v2 = - u2 / r
          u2 = v2 / v1
!
          DO 110 j = k, n
            t = a (l, j) + u2 * a (l1, j)
            a (l, j) = a (l, j) + t * v1
            a (l1, j) = a (l1, j) + t * v2
  110     END DO
!
          a (l1, k) = 0.0d0
!
          DO 120 j = l, n
            t = b (l, j) + u2 * b (l1, j)
            b (l, j) = b (l, j) + t * v1
            b (l1, j) = b (l1, j) + t * v2
  120     END DO
!     .......... zero b(l+1,l) ..........
          s = dabs (b (l1, l1) ) + dabs (b (l1, l) )
          IF (s.eq.0.0d0) goto 150
          u1 = b (l1, l1) / s
          u2 = b (l1, l) / s
          r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
          v1 = - (u1 + r) / r
          v2 = - u2 / r
          u2 = v2 / v1
!
          DO 130 i = 1, l1
            t = b (i, l1) + u2 * b (i, l)
            b (i, l1) = b (i, l1) + t * v1
            b (i, l) = b (i, l) + t * v2
  130     END DO
!
          b (l1, l) = 0.0d0
!
          DO 140 i = 1, n
            t = a (i, l1) + u2 * a (i, l)
            a (i, l1) = a (i, l1) + t * v1
            a (i, l) = a (i, l) + t * v2
  140     END DO
!
          IF (.not.matz) goto 150
!
          DO 145 i = 1, n
            t = z (i, l1) + u2 * z (i, l)
            z (i, l1) = z (i, l1) + t * v1
            z (i, l) = z (i, l) + t * v2
  145     END DO
!
  150   END DO
!
  160 END DO
!
  170 RETURN
      END SUBROUTINE qzhes

      SUBROUTINE qzit (nm, n, a, b, eps1, matz, z, ierr)
!
      INTEGER i, j, k, l, n, en, k1, k2, ld, ll, l1, na, nm, ish, itn,  &
      its, km1, lm1, enm2, ierr, lor1, enorn
      DOUBLE PRECISION a (nm, n), b (nm, n), z (nm, n)
      DOUBLE PRECISION r, s, t, a1, a2, a3, ep, sh, u1, u2, u3, v1, v2,  &
      v3, ani, a11, a12, a21, a22, a33, a34, a43, a44, bni, b11, b12,   &
      b22, b33, b34, b44, epsa, epsb, eps1, anorm, bnorm, epslon
      LOGICAL matz, notlas
!
!     this subroutine is the second step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
!     as modified in technical note nasa tn d-7305(1973) by ward.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in upper hessenberg form and the other in upper triangular form.
!     it reduces the hessenberg matrix to quasi-triangular form using
!     orthogonal transformations while maintaining the triangular form
!     of the other matrix.  it is usually preceded by  qzhes  and
!     followed by  qzval  and, possibly,  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper hessenberg matrix.
!
!        b contains a real upper triangular matrix.
!
!        eps1 is a tolerance used to determine negligible elements.
!          eps1 = 0.0 (or negative) may be input, in which case an
!          element will be neglected only if it is less than roundoff
!          error times the norm of its matrix.  if the input eps1 is
!          positive, then an element will be considered negligible
!          if it is less than eps1 times the norm of its matrix.  a
!          positive value of eps1 may result in faster execution,
!          but less accurate results.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reduction
!          by  qzhes, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced to quasi-triangular form.  the elements
!          below the first subdiagonal are still zero and no two
!          consecutive subdiagonal elements are nonzero.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  the location b(n,1) is used to store
!          eps1 times the norm of b for later use by  qzval  and  qzvec.
!
!        z contains the product of the right hand transformations
!          (for both steps) if matz has been set to .true..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
!     .......... compute epsa,epsb ..........
      anorm = 0.0d0
      bnorm = 0.0d0
!
      DO 30 i = 1, n
        ani = 0.0d0
        IF (i.ne.1) ani = dabs (a (i, i - 1) )
        bni = 0.0d0
!
        DO 20 j = i, n
          ani = ani + dabs (a (i, j) )
          bni = bni + dabs (b (i, j) )
   20   END DO
!
        IF (ani.gt.anorm) anorm = ani
        IF (bni.gt.bnorm) bnorm = bni
   30 END DO
!
      IF (anorm.eq.0.0d0) anorm = 1.0d0
      IF (bnorm.eq.0.0d0) bnorm = 1.0d0
      ep = eps1
      IF (ep.gt.0.0d0) goto 50
!     .......... use roundoff level if eps1 is zero ..........
      epslon = 1.0e0
      ep = epslon
   50 epsa = ep * anorm
      epsb = ep * bnorm
!     .......... reduce a to quasi-triangular form, while
!                keeping b triangular ..........
      lor1 = 1
      enorn = n
      en = n
      itn = 30 * n
!     .......... begin qz step ..........
   60 IF (en.le.2) goto 1001
      IF (.not.matz) enorn = en
      its = 0
      na = en - 1
      enm2 = na - 1
   70 ish = 2
!     .......... check for convergence or reducibility.
!                for l=en step -1 until 1 do -- ..........
      DO 80 ll = 1, en
        lm1 = en - ll
        l = lm1 + 1
        IF (l.eq.1) goto 95
        IF (dabs (a (l, lm1) ) .le.epsa) goto 90
   80 END DO
!
   90 a (l, lm1) = 0.0d0
      IF (l.lt.na) goto 95
!     .......... 1-by-1 or 2-by-2 block isolated ..........
      en = lm1
      GOTO 60
!     .......... check for small top of b ..........
   95 ld = l
  100 l1 = l + 1
      b11 = b (l, l)
      IF (dabs (b11) .gt.epsb) goto 120
      b (l, l) = 0.0d0
      s = dabs (a (l, l) ) + dabs (a (l1, l) )
      u1 = a (l, l) / s
      u2 = a (l1, l) / s
      r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
      v1 = - (u1 + r) / r
      v2 = - u2 / r
      u2 = v2 / v1
!
      DO 110 j = l, enorn
        t = a (l, j) + u2 * a (l1, j)
        a (l, j) = a (l, j) + t * v1
        a (l1, j) = a (l1, j) + t * v2
        t = b (l, j) + u2 * b (l1, j)
        b (l, j) = b (l, j) + t * v1
        b (l1, j) = b (l1, j) + t * v2
  110 END DO
!
      IF (l.ne.1) a (l, lm1) = - a (l, lm1)
      lm1 = l
      l = l1
      GOTO 90
  120 a11 = a (l, l) / b11
      a21 = a (l1, l) / b11
      IF (ish.eq.1) goto 140
!     .......... iteration strategy ..........
      IF (itn.eq.0) goto 1000
      IF (its.eq.10) goto 155
!     .......... determine type of shift ..........
      b22 = b (l1, l1)
      IF (dabs (b22) .lt.epsb) b22 = epsb
      b33 = b (na, na)
      IF (dabs (b33) .lt.epsb) b33 = epsb
      b44 = b (en, en)
      IF (dabs (b44) .lt.epsb) b44 = epsb
      a33 = a (na, na) / b33
      a34 = a (na, en) / b44
      a43 = a (en, na) / b33
      a44 = a (en, en) / b44
      b34 = b (na, en) / b44
      t = 0.5d0 * (a43 * b34 - a33 - a44)
      r = t * t + a34 * a43 - a33 * a44
      IF (r.lt.0.0d0) goto 150
!     .......... determine single shift zeroth column of a ..........
      ish = 1
      r = dsqrt (r)
      sh = - t + r
      s = - t - r
      IF (dabs (s - a44) .lt.dabs (sh - a44) ) sh = s
!     .......... look for two consecutive small
!                sub-diagonal elements of a.
!                for l=en-2 step -1 until ld do -- ..........
      DO 130 ll = ld, enm2
        l = enm2 + ld-ll
        IF (l.eq.ld) goto 140
        lm1 = l - 1
        l1 = l + 1
        t = a (l, l)
        IF (dabs (b (l, l) ) .gt.epsb) t = t - sh * b (l, l)
        IF (dabs (a (l, lm1) ) .le.dabs (t / a (l1, l) ) * epsa) goto   &
        100
  130 END DO
!
  140 a1 = a11 - sh
      a2 = a21
      IF (l.ne.ld) a (l, lm1) = - a (l, lm1)
      GOTO 160
!     .......... determine double shift zeroth column of a ..........
  150 a12 = a (l, l1) / b22
      a22 = a (l1, l1) / b22
      b12 = b (l, l1) / b22
      a1 = ( (a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)   &
      / a21 + a12 - a11 * b12
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 *  &
      b34
      a3 = a (l1 + 1, l1) / b22
      GOTO 160
!     .......... ad hoc shift ..........
  155 a1 = 0.0d0
      a2 = 1.0d0
      a3 = 1.1605d0
  160 its = its + 1
      itn = itn - 1
      IF (.not.matz) lor1 = ld
!     .......... main loop ..........
      DO 260 k = l, na
        notlas = k.ne.na.and.ish.eq.2
        k1 = k + 1
        k2 = k + 2
        km1 = max0 (k - 1, l)
        ll = min0 (en, k1 + ish)
        IF (notlas) goto 190
!     .......... zero a(k+1,k-1) ..........
        IF (k.eq.l) goto 170
        a1 = a (k, km1)
        a2 = a (k1, km1)
  170   s = dabs (a1) + dabs (a2)
        IF (s.eq.0.0d0) goto 70
        u1 = a1 / s
        u2 = a2 / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        u2 = v2 / v1
!
        DO 180 j = km1, enorn
          t = a (k, j) + u2 * a (k1, j)
          a (k, j) = a (k, j) + t * v1
          a (k1, j) = a (k1, j) + t * v2
          t = b (k, j) + u2 * b (k1, j)
          b (k, j) = b (k, j) + t * v1
          b (k1, j) = b (k1, j) + t * v2
  180   END DO
!
        IF (k.ne.l) a (k1, km1) = 0.0d0
        GOTO 240
!     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
  190   IF (k.eq.l) goto 200
        a1 = a (k, km1)
        a2 = a (k1, km1)
        a3 = a (k2, km1)
  200   s = dabs (a1) + dabs (a2) + dabs (a3)
        IF (s.eq.0.0d0) goto 260
        u1 = a1 / s
        u2 = a2 / s
        u3 = a3 / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2 + u3 * u3), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        v3 = - u3 / r
        u2 = v2 / v1
        u3 = v3 / v1
!
        DO 210 j = km1, enorn
          t = a (k, j) + u2 * a (k1, j) + u3 * a (k2, j)
          a (k, j) = a (k, j) + t * v1
          a (k1, j) = a (k1, j) + t * v2
          a (k2, j) = a (k2, j) + t * v3
          t = b (k, j) + u2 * b (k1, j) + u3 * b (k2, j)
          b (k, j) = b (k, j) + t * v1
          b (k1, j) = b (k1, j) + t * v2
          b (k2, j) = b (k2, j) + t * v3
  210   END DO
!
        IF (k.eq.l) goto 220
        a (k1, km1) = 0.0d0
        a (k2, km1) = 0.0d0
!     .......... zero b(k+2,k+1) and b(k+2,k) ..........
  220   s = dabs (b (k2, k2) ) + dabs (b (k2, k1) ) + dabs (b (k2, k) )
        IF (s.eq.0.0d0) goto 240
        u1 = b (k2, k2) / s
        u2 = b (k2, k1) / s
        u3 = b (k2, k) / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2 + u3 * u3), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        v3 = - u3 / r
        u2 = v2 / v1
        u3 = v3 / v1
!
        DO 230 i = lor1, ll
          t = a (i, k2) + u2 * a (i, k1) + u3 * a (i, k)
          a (i, k2) = a (i, k2) + t * v1
          a (i, k1) = a (i, k1) + t * v2
          a (i, k) = a (i, k) + t * v3
          t = b (i, k2) + u2 * b (i, k1) + u3 * b (i, k)
          b (i, k2) = b (i, k2) + t * v1
          b (i, k1) = b (i, k1) + t * v2
          b (i, k) = b (i, k) + t * v3
  230   END DO
!
        b (k2, k) = 0.0d0
        b (k2, k1) = 0.0d0
        IF (.not.matz) goto 240
!
        DO 235 i = 1, n
          t = z (i, k2) + u2 * z (i, k1) + u3 * z (i, k)
          z (i, k2) = z (i, k2) + t * v1
          z (i, k1) = z (i, k1) + t * v2
          z (i, k) = z (i, k) + t * v3
  235   END DO
!     .......... zero b(k+1,k) ..........
  240   s = dabs (b (k1, k1) ) + dabs (b (k1, k) )
        IF (s.eq.0.0d0) goto 260
        u1 = b (k1, k1) / s
        u2 = b (k1, k) / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        u2 = v2 / v1
!
        DO 250 i = lor1, ll
          t = a (i, k1) + u2 * a (i, k)
          a (i, k1) = a (i, k1) + t * v1
          a (i, k) = a (i, k) + t * v2
          t = b (i, k1) + u2 * b (i, k)
          b (i, k1) = b (i, k1) + t * v1
          b (i, k) = b (i, k) + t * v2
  250   END DO
!
        b (k1, k) = 0.0d0
        IF (.not.matz) goto 260
!
        DO 255 i = 1, n
          t = z (i, k1) + u2 * z (i, k)
          z (i, k1) = z (i, k1) + t * v1
          z (i, k) = z (i, k) + t * v2
  255   END DO
!
  260 END DO
!     .......... end qz step ..........
      GOTO 70
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
!     .......... save epsb for use by qzval and qzvec ..........
 1001 IF (n.gt.1) b (n, 1) = epsb
      RETURN
      END SUBROUTINE qzit

      SUBROUTINE qzval (nm, n, a, b, alfr, alfi, beta, matz, z)
!
      INTEGER i, j, n, en, na, nm, nn, isw
      DOUBLE PRECISION a (nm, n), b (nm, n), alfr (n), alfi (n), beta (n)&
      , z (nm, n)
      DOUBLE PRECISION c, d, e, r, s, t, an, a1, a2, bn, cq, cz, di, dr, &
      ei, ti, tr, u1, u2, v1, v2, a1i, a11, a12, a2i, a21, a22, b11,    &
      b12, b22, sqi, sqr, ssi, ssr, szi, szr, a11i, a11r, a12i, a12r,   &
      a22i, a22r, epsb
      LOGICAL matz
!
!     this subroutine is the third step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them
!     in quasi-triangular form and the other in upper triangular form.
!     it reduces the quasi-triangular matrix further, so that any
!     remaining 2-by-2 blocks correspond to pairs of complex
!     eigenvalues, and returns quantities whose ratios give the
!     generalized eigenvalues.  it is usually preceded by  qzhes
!     and  qzit  and may be followed by  qzvec.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        matz should be set to .true. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .false. otherwise.
!
!        z contains, if matz has been set to .true., the
!          transformation matrix produced in the reductions by qzhes
!          and qzit, if performed, or else the identity matrix.
!          if matz has been set to .false., z is not referenced.
!
!     on output
!
!        a has been reduced further to a quasi-triangular matrix
!          in which all nonzero subdiagonal elements correspond to
!          pairs of complex eigenvalues.
!
!        b is still in upper triangular form, although its elements
!          have been altered.  b(n,1) is unaltered.
!
!        alfr and alfi contain the real and imaginary parts of the
!          diagonal elements of the triangular matrix that would be
!          obtained if a were reduced completely to triangular form
!          by unitary transformations.  non-zero values of alfi occur
!          in pairs, the first member positive and the second negative.
!
!        beta contains the diagonal elements of the corresponding b,
!          normalized to be real and non-negative.  the generalized
!          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
!
!        z contains the product of the right hand transformations
!          (for all three steps) if matz has been set to .true.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      epsb = b (n, 1)
      isw = 1
!     .......... find eigenvalues of quasi-triangular matrices.
!                for en=n step -1 until 1 do -- ..........
      DO 510 nn = 1, n
        en = n + 1 - nn
        na = en - 1
        IF (isw.eq.2) goto 505
        IF (en.eq.1) goto 410
        IF (a (en, na) .ne.0.0d0) goto 420
!     .......... 1-by-1 block, one real root ..........
  410   alfr (en) = a (en, en)
        IF (b (en, en) .lt.0.0d0) alfr (en) = - alfr (en)
        beta (en) = dabs (b (en, en) )
        alfi (en) = 0.0d0
        GOTO 510
!     .......... 2-by-2 block ..........
  420   IF (dabs (b (na, na) ) .le.epsb) goto 455
        IF (dabs (b (en, en) ) .gt.epsb) goto 430
        a1 = a (en, en)
        a2 = a (en, na)
        bn = 0.0d0
        GOTO 435
  430   an = dabs (a (na, na) ) + dabs (a (na, en) ) + dabs (a (en, na) &
        ) + dabs (a (en, en) )
        bn = dabs (b (na, na) ) + dabs (b (na, en) ) + dabs (b (en, en) &
        )
        a11 = a (na, na) / an
        a12 = a (na, en) / an
        a21 = a (en, na) / an
        a22 = a (en, en) / an
        b11 = b (na, na) / bn
        b12 = b (na, en) / bn
        b22 = b (en, en) / bn
        e = a11 / b11
        ei = a22 / b22
        s = a21 / (b11 * b22)
        t = (a22 - e * b22) / b22
        IF (dabs (e) .le.dabs (ei) ) goto 431
        e = ei
        t = (a11 - e * b11) / b11
  431   c = 0.5d0 * (t - s * b12)
        d = c * c + s * (a12 - e * b12)
        IF (d.lt.0.0d0) goto 480
!     .......... two real roots.
!                zero both a(en,na) and b(en,na) ..........
        e = e+ (c + dsign (dsqrt (d), c) )
        a11 = a11 - e * b11
        a12 = a12 - e * b12
        a22 = a22 - e * b22
        IF (dabs (a11) + dabs (a12) .lt.dabs (a21) + dabs (a22) ) goto  &
        432
        a1 = a12
        a2 = a11
        GOTO 435
  432   a1 = a22
        a2 = a21
!     .......... choose and apply real z ..........
  435   s = dabs (a1) + dabs (a2)
        u1 = a1 / s
        u2 = a2 / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        u2 = v2 / v1
!
        DO 440 i = 1, en
          t = a (i, en) + u2 * a (i, na)
          a (i, en) = a (i, en) + t * v1
          a (i, na) = a (i, na) + t * v2
          t = b (i, en) + u2 * b (i, na)
          b (i, en) = b (i, en) + t * v1
          b (i, na) = b (i, na) + t * v2
  440   END DO
!
        IF (.not.matz) goto 450
!
        DO 445 i = 1, n
          t = z (i, en) + u2 * z (i, na)
          z (i, en) = z (i, en) + t * v1
          z (i, na) = z (i, na) + t * v2
  445   END DO
!
  450   IF (bn.eq.0.0d0) goto 475
        IF (an.lt.dabs (e) * bn) goto 455
        a1 = b (na, na)
        a2 = b (en, na)
        GOTO 460
  455   a1 = a (na, na)
        a2 = a (en, na)
!     .......... choose and apply real q ..........
  460   s = dabs (a1) + dabs (a2)
        IF (s.eq.0.0d0) goto 475
        u1 = a1 / s
        u2 = a2 / s
        r = dsign (dsqrt (u1 * u1 + u2 * u2), u1)
        v1 = - (u1 + r) / r
        v2 = - u2 / r
        u2 = v2 / v1
!
        DO 470 j = na, n
          t = a (na, j) + u2 * a (en, j)
          a (na, j) = a (na, j) + t * v1
          a (en, j) = a (en, j) + t * v2
          t = b (na, j) + u2 * b (en, j)
          b (na, j) = b (na, j) + t * v1
          b (en, j) = b (en, j) + t * v2
  470   END DO
!
  475   a (en, na) = 0.0d0
        b (en, na) = 0.0d0
        alfr (na) = a (na, na)
        alfr (en) = a (en, en)
        IF (b (na, na) .lt.0.0d0) alfr (na) = - alfr (na)
        IF (b (en, en) .lt.0.0d0) alfr (en) = - alfr (en)
        beta (na) = dabs (b (na, na) )
        beta (en) = dabs (b (en, en) )
        alfi (en) = 0.0d0
        alfi (na) = 0.0d0
        GOTO 505
!     .......... two complex roots ..........
  480   e = e+c
        ei = dsqrt ( - d)
        a11r = a11 - e * b11
        a11i = ei * b11
        a12r = a12 - e * b12
        a12i = ei * b12
        a22r = a22 - e * b22
        a22i = ei * b22
        IF (dabs (a11r) + dabs (a11i) + dabs (a12r) + dabs (a12i)       &
        .lt.dabs (a21) + dabs (a22r) + dabs (a22i) ) goto 482
        a1 = a12r
        a1i = a12i
        a2 = - a11r
        a2i = - a11i
        GOTO 485
  482   a1 = a22r
        a1i = a22i
        a2 = - a21
        a2i = 0.0d0
!     .......... choose complex z ..........
  485   cz = dsqrt (a1 * a1 + a1i * a1i)
        IF (cz.eq.0.0d0) goto 487
        szr = (a1 * a2 + a1i * a2i) / cz
        szi = (a1 * a2i - a1i * a2) / cz
        r = dsqrt (cz * cz + szr * szr + szi * szi)
        cz = cz / r
        szr = szr / r
        szi = szi / r
        GOTO 490
  487   szr = 1.0d0
        szi = 0.0d0
  490   IF (an.lt. (dabs (e) + ei) * bn) goto 492
        a1 = cz * b11 + szr * b12
        a1i = szi * b12
        a2 = szr * b22
        a2i = szi * b22
        GOTO 495
  492   a1 = cz * a11 + szr * a12
        a1i = szi * a12
        a2 = cz * a21 + szr * a22
        a2i = szi * a22
!     .......... choose complex q ..........
  495   cq = dsqrt (a1 * a1 + a1i * a1i)
        IF (cq.eq.0.0d0) goto 497
        sqr = (a1 * a2 + a1i * a2i) / cq
        sqi = (a1 * a2i - a1i * a2) / cq
        r = dsqrt (cq * cq + sqr * sqr + sqi * sqi)
        cq = cq / r
        sqr = sqr / r
        sqi = sqi / r
        GOTO 500
  497   sqr = 1.0d0
        sqi = 0.0d0
!     .......... compute diagonal elements that would result
!                if transformations were applied ..........
  500   ssr = sqr * szr + sqi * szi
        ssi = sqr * szi - sqi * szr
        i = 1
        tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr *    &
        a22
        ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
        dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
        di = cq * szi * b12 + ssi * b22
        GOTO 503
  502   i = 2
        tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz *    &
        a22
        ti = - ssi * a11 - sqi * cz * a12 + cq * szi * a21
        dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
        di = - ssi * b11 - sqi * cz * b12
  503   t = ti * dr - tr * di
        j = na
        IF (t.lt.0.0d0) j = en
        r = dsqrt (dr * dr + di * di)
        beta (j) = bn * r
        alfr (j) = an * (tr * dr + ti * di) / r
        alfi (j) = an * t / r
        IF (i.eq.1) goto 502
  505   isw = 3 - isw
  510 END DO
      b (n, 1) = epsb
!
      RETURN
      END SUBROUTINE qzval

      SUBROUTINE qzvec (nm, n, a, b, alfr, alfi, beta, z)
!
      INTEGER i, j, k, m, n, en, ii, jj, na, nm, nn, isw, enm2
      DOUBLE PRECISION a (nm, n), b (nm, n), alfr (n), alfi (n), beta (n)&
      , z (nm, n)
      DOUBLE PRECISION d, q, r, s, t, w, x, y, di, dr, ra, rr, sa, ti,   &
      tr, t1, t2, w1, x1, zz, z1, alfm, almi, almr, betm, epsb
!
!     this subroutine is the optional fourth step of the qz algorithm
!     for solving generalized matrix eigenvalue problems,
!     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
!
!     this subroutine accepts a pair of real matrices, one of them in
!     quasi-triangular form (in which each 2-by-2 block corresponds to
!     a pair of complex eigenvalues) and the other in upper triangular
!     form.  it computes the eigenvectors of the triangular problem and
!     transforms the results back to the original coordinate system.
!     it is usually preceded by  qzhes,  qzit, and  qzval.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices.
!
!        a contains a real upper quasi-triangular matrix.
!
!        b contains a real upper triangular matrix.  in addition,
!          location b(n,1) contains the tolerance quantity (epsb)
!          computed and saved in  qzit.
!
!        alfr, alfi, and beta  are vectors with components whose
!          ratios ((alfr+i*alfi)/beta) are the generalized
!          eigenvalues.  they are usually obtained from  qzval.
!
!        z contains the transformation matrix produced in the
!          reductions by  qzhes,  qzit, and  qzval, if performed.
!          if the eigenvectors of the triangular problem are
!          desired, z must contain the identity matrix.
!
!     on output
!
!        a is unaltered.  its subdiagonal elements provide information
!           about the storage of the complex eigenvectors.
!
!        b has been destroyed.
!
!        alfr, alfi, and beta are unaltered.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
!            the i-th column of z contains its eigenvector.
!          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
!            if alfi(i) .gt. 0.0, the eigenvalue is the first of
!              a complex pair and the i-th and (i+1)-th columns
!              of z contain its eigenvector.
!            if alfi(i) .lt. 0.0, the eigenvalue is the second of
!              a complex pair and the (i-1)-th and i-th columns
!              of z contain the conjugate of its eigenvector.
!          each eigenvector is normalized so that the modulus
!          of its largest component is 1.0 .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      epsb = b (n, 1)
      isw = 1
!     .......... for en=n step -1 until 1 do -- ..........
      DO 800 nn = 1, n
        en = n + 1 - nn
        na = en - 1
        IF (isw.eq.2) goto 795
        IF (alfi (en) .ne.0.0d0) goto 710
!     .......... real vector ..........
        m = en
        b (en, en) = 1.0d0
        IF (na.eq.0) goto 800
        alfm = alfr (m)
        betm = beta (m)
!     .......... for i=en-1 step -1 until 1 do -- ..........
        DO 700 ii = 1, na
          i = en - ii
          w = betm * a (i, i) - alfm * b (i, i)
          r = 0.0d0
!
          DO 610 j = m, en
  610     r = r + (betm * a (i, j) - alfm * b (i, j) ) * b (j, en)
!
          IF (i.eq.1.or.isw.eq.2) goto 630
          IF (betm * a (i, i - 1) .eq.0.0d0) goto 630
          zz = w
          s = r
          GOTO 690
  630     m = i
          IF (isw.eq.2) goto 640
!     .......... real 1-by-1 block ..........
          t = w
          IF (w.eq.0.0d0) t = epsb
          b (i, en) = - r / t
          GOTO 700
!     .......... real 2-by-2 block ..........
  640     x = betm * a (i, i + 1) - alfm * b (i, i + 1)
          y = betm * a (i + 1, i)
          q = w * zz - x * y
          t = (x * s - zz * r) / q
          b (i, en) = t
          IF (dabs (x) .le.dabs (zz) ) goto 650
          b (i + 1, en) = ( - r - w * t) / x
          GOTO 690
  650     b (i + 1, en) = ( - s - y * t) / zz
  690     isw = 3 - isw
  700   END DO
!     .......... end real vector ..........
        GOTO 800
!     .......... complex vector ..........
  710   m = na
        almr = alfr (m)
        almi = alfi (m)
        betm = beta (m)
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
        y = betm * a (en, na)
        b (na, na) = - almi * b (en, en) / y
        b (na, en) = (almr * b (en, en) - betm * a (en, en) ) / y
        b (en, na) = 0.0d0
        b (en, en) = 1.0d0
        enm2 = na - 1
        IF (enm2.eq.0) goto 795
!     .......... for i=en-2 step -1 until 1 do -- ..........
        DO 790 ii = 1, enm2
          i = na - ii
          w = betm * a (i, i) - almr * b (i, i)
          w1 = - almi * b (i, i)
          ra = 0.0d0
          sa = 0.0d0
!
          DO 760 j = m, en
            x = betm * a (i, j) - almr * b (i, j)
            x1 = - almi * b (i, j)
            ra = ra + x * b (j, na) - x1 * b (j, en)
            sa = sa + x * b (j, en) + x1 * b (j, na)
  760     END DO
!
          IF (i.eq.1.or.isw.eq.2) goto 770
          IF (betm * a (i, i - 1) .eq.0.0d0) goto 770
          zz = w
          z1 = w1
          r = ra
          s = sa
          isw = 2
          GOTO 790
  770     m = i
          IF (isw.eq.2) goto 780
!     .......... complex 1-by-1 block ..........
          tr = - ra
          ti = - sa
  773     dr = w
          di = w1
!     .......... complex divide (t1,t2) = (tr,ti) / (dr,di) ..........
  775     IF (dabs (di) .gt.dabs (dr) ) goto 777
          rr = di / dr
          d = dr + di * rr
          t1 = (tr + ti * rr) / d
          t2 = (ti - tr * rr) / d
          GOTO (787, 782), isw
  777     rr = dr / di
          d = dr * rr + di
          t1 = (tr * rr + ti) / d
          t2 = (ti * rr - tr) / d
          GOTO (787, 782), isw
!     .......... complex 2-by-2 block ..........
  780     x = betm * a (i, i + 1) - almr * b (i, i + 1)
          x1 = - almi * b (i, i + 1)
          y = betm * a (i + 1, i)
          tr = y * ra - w * r + w1 * s
          ti = y * sa - w * s - w1 * r
          dr = w * zz - w1 * z1 - x * y
          di = w * z1 + w1 * zz - x1 * y
          IF (dr.eq.0.0d0.and.di.eq.0.0d0) dr = epsb
          GOTO 775
  782     b (i + 1, na) = t1
          b (i + 1, en) = t2
          isw = 1
          IF (dabs (y) .gt.dabs (w) + dabs (w1) ) goto 785
          tr = - ra - x * b (i + 1, na) + x1 * b (i + 1, en)
          ti = - sa - x * b (i + 1, en) - x1 * b (i + 1, na)
          GOTO 773
  785     t1 = ( - r - zz * b (i + 1, na) + z1 * b (i + 1, en) )        &
          / y
          t2 = ( - s - zz * b (i + 1, en) - z1 * b (i + 1, na) )        &
          / y
  787     b (i, na) = t1
          b (i, en) = t2
  790   END DO
!     .......... end complex vector ..........
  795   isw = 3 - isw
  800 END DO
!     .......... end back substitution.
!                transform to original coordinate system.
!                for j=n step -1 until 1 do -- ..........
      DO 880 jj = 1, n
        j = n + 1 - jj
!
        DO 880 i = 1, n
          zz = 0.0d0
!
          DO 860 k = 1, j
  860     zz = zz + z (i, k) * b (k, j)
!
          z (i, j) = zz
  880 CONTINUE
!     .......... normalize so that modulus of largest
!                component of each vector is 1.
!                (isw is 1 initially from before) ..........
      DO 950 j = 1, n
        d = 0.0d0
        IF (isw.eq.2) goto 920
        IF (alfi (j) .ne.0.0d0) goto 945
!
        DO 890 i = 1, n
          IF (dabs (z (i, j) ) .gt.d) d = dabs (z (i, j) )
  890   END DO
!
        DO 900 i = 1, n
  900   z (i, j) = z (i, j) / d
!
        GOTO 950
!
  920   DO 930 i = 1, n
          r = dabs (z (i, j - 1) ) + dabs (z (i, j) )
          IF (r.ne.0.0d0) r = r * dsqrt ( (z (i, j - 1) / r) **2 +      &
          (z (i, j) / r) **2)
          IF (r.gt.d) d = r
  930   END DO
!
        DO 940 i = 1, n
          z (i, j - 1) = z (i, j - 1) / d
          z (i, j) = z (i, j) / d
  940   END DO
!
  945   isw = 3 - isw
  950 END DO
!
      RETURN
      END SUBROUTINE qzvec

      SUBROUTINE rgg (nm, n, a, b, alfr, alfi, beta, matz, z, ierr)
!
      INTEGER n, nm, ierr, matz
      DOUBLE PRECISION a(nm, n), b(nm, n), alfr(n), alfi(n), beta(n), z(nm, n)
      LOGICAL tf
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real general generalized eigenproblem  ax = (lambda)bx.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real general matrix.
!
!        b  contains a real general matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        alfr  and  alfi  contain the real and imaginary parts,
!        respectively, of the numerators of the eigenvalues.
!
!        beta  contains the denominators of the eigenvalues,
!        which are thus given by the ratios  (alfr+i*alfi)/beta.
!        complex conjugate pairs of eigenvalues appear consecutively
!        with the eigenvalue having the positive imaginary part first.
!
!        z  contains the real and imaginary parts of the eigenvectors
!        if matz is not zero.  if the j-th eigenvalue is real, the
!        j-th column of  z  contains its eigenvector.  if the j-th
!        eigenvalue is complex with positive imaginary part, the
!        j-th and (j+1)-th columns of  z  contain the real and
!        imaginary parts of its eigenvector.  the conjugate of this
!        vector is the eigenvector for the conjugate eigenvalue.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for qzit.
!           the normal completion code is zero.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!

      IF (n.le.nm) goto 10
      ierr = 10 * n
      GOTO 50
!
   10 IF (matz.ne.0) goto 20
!     .......... find eigenvalues only ..........
      tf = .false.
      CALL qzhes (nm, n, a, b, tf, z)
      CALL qzit (nm, n, a, b, 1.0e-14_8, tf, z, ierr)
      CALL qzval (nm, n, a, b, alfr, alfi, beta, tf, z)
      GOTO 50
!     .......... find both eigenvalues and eigenvectors ..........
   20 tf = .true.
      CALL qzhes (nm, n, a, b, tf, z)
      CALL qzit (nm, n, a, b, 1.0e-14_8, tf, z, ierr)
      CALL qzval (nm, n, a, b, alfr, alfi, beta, tf, z)
      IF (ierr.ne.0) goto 50
      CALL qzvec (nm, n, a, b, alfr, alfi, beta, z)
   50 RETURN
      END SUBROUTINE rgg

      SUBROUTINE hqr2 (nm, n, low, igh, h, wr, wi, z, ierr)
!
      INTEGER i, j, k, l, m, n, en, ii, jj, ll, mm, na, nm, nn, igh,    &
      itn, its, low, mp2, enm2, ierr
      DOUBLEPRECISION h (nm, n), wr (n), wi (n), z (nm, n)
      DOUBLEPRECISION p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, norm, &
      tst1, tst2
      LOGICAL notlas
!
!     this subroutine is a translation of the algol procedure hqr2,
!     num. math. 16, 181-204(1970) by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a real upper hessenberg matrix by the qr method.  the
!     eigenvectors of a real general matrix can also be found
!     if  elmhes  and  eltran  or  orthes  and  ortran  have
!     been used to reduce this general matrix to hessenberg form
!     and to accumulate the similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        low and igh are integers determined by the balancing
!          subroutine  balanc.  if  balanc  has not been used,
!          set low=1, igh=n.
!
!        h contains the upper hessenberg matrix.
!
!        z contains the transformation matrix produced by  eltran
!          after the reduction by  elmhes, or by  ortran  after the
!          reduction by  orthes, if performed.  if the eigenvectors
!          of the hessenberg matrix are desired, z must contain the
!          identity matrix.
!
!     on output
!
!        h has been destroyed.
!
!        wr and wi contain the real and imaginary parts,
!          respectively, of the eigenvalues.  the eigenvalues
!          are unordered except that complex conjugate pairs
!          of values appear consecutively with the eigenvalue
!          having the positive imaginary part first.  if an
!          error exit is made, the eigenvalues should be correct
!          for indices ierr+1,...,n.
!
!        z contains the real and imaginary parts of the eigenvectors.
!          if the i-th eigenvalue is real, the i-th column of z
!          contains its eigenvector.  if the i-th eigenvalue is complex
!          with positive imaginary part, the i-th and (i+1)-th
!          columns of z contain the real and imaginary parts of its
!          eigenvector.  the eigenvectors are unnormalized.  if an
!          error exit is made, none of the eigenvectors has been found.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the limit of 30*n iterations is exhausted
!                     while the j-th eigenvalue is being sought.
!
!     calls cdiv for complex division.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
      ierr = 0
      norm = 0.0d0
      k = 1
!     .......... store roots isolated by balanc
!                and compute matrix norm ..........
      DO 50 i = 1, n
!
        DO 40 j = k, n
   40   norm = norm + dabs (h (i, j) )
!
        k = i
        IF (i.ge.low.and.i.le.igh) goto 50
        wr (i) = h (i, i)
        wi (i) = 0.0d0
   50 END DO
!
      en = igh
      t = 0.0d0
      itn = 30 * n
!     .......... search for next eigenvalues ..........
   60 IF (en.lt.low) goto 340
      its = 0
      na = en - 1
      enm2 = na - 1
!     .......... look for single small sub-diagonal element
!                for l=en step -1 until low do -- ..........
   70 DO 80 ll = low, en
        l = en + low - ll
        IF (l.eq.low) goto 100
        s = dabs (h (l - 1, l - 1) ) + dabs (h (l, l) )
        IF (s.eq.0.0d0) s = norm
        tst1 = s
        tst2 = tst1 + dabs (h (l, l - 1) )
        IF (tst2.eq.tst1) goto 100
   80 END DO
!     .......... form shift ..........
  100 x = h (en, en)
      IF (l.eq.en) goto 270
      y = h (na, na)
      w = h (en, na) * h (na, en)
      IF (l.eq.na) goto 280
      IF (itn.eq.0) goto 1000
      IF (its.ne.10.and.its.ne.20) goto 130
!     .......... form exceptional shift ..........
      t = t + x
!
      DO 120 i = low, en
  120 h (i, i) = h (i, i) - x
!
      s = dabs (h (en, na) ) + dabs (h (na, enm2) )
      x = 0.75d0 * s
      y = x
      w = - 0.4375d0 * s * s
  130 its = its + 1
      itn = itn - 1
!     .......... look for two consecutive small
!                sub-diagonal elements.
!                for m=en-2 step -1 until l do -- ..........
      DO 140 mm = l, enm2
        m = enm2 + l - mm
        zz = h (m, m)
        r = x - zz
        s = y - zz
        p = (r * s - w) / h (m + 1, m) + h (m, m + 1)
        q = h (m + 1, m + 1) - zz - r - s
        r = h (m + 2, m + 1)
        s = dabs (p) + dabs (q) + dabs (r)
        p = p / s
        q = q / s
        r = r / s
        IF (m.eq.l) goto 150
        tst1 = dabs (p) * (dabs (h (m - 1, m - 1) ) + dabs (zz) + dabs (&
        h (m + 1, m + 1) ) )
        tst2 = tst1 + dabs (h (m, m - 1) ) * (dabs (q) + dabs (r) )
        IF (tst2.eq.tst1) goto 150
  140 END DO
!
  150 mp2 = m + 2
!
      DO 160 i = mp2, en
        h (i, i - 2) = 0.0d0
        IF (i.eq.mp2) goto 160
        h (i, i - 3) = 0.0d0
  160 END DO
!     .......... double qr step involving rows l to en and
!                columns m to en ..........
      DO 260 k = m, na
        notlas = k.ne.na
        IF (k.eq.m) goto 170
        p = h (k, k - 1)
        q = h (k + 1, k - 1)
        r = 0.0d0
        IF (notlas) r = h (k + 2, k - 1)
        x = dabs (p) + dabs (q) + dabs (r)
        IF (x.eq.0.0d0) goto 260
        p = p / x
        q = q / x
        r = r / x
  170   s = dsign (dsqrt (p * p + q * q + r * r), p)
        IF (k.eq.m) goto 180
        h (k, k - 1) = - s * x
        GOTO 190
  180   IF (l.ne.m) h (k, k - 1) = - h (k, k - 1)
  190   p = p + s
        x = p / s
        y = q / s
        zz = r / s
        q = q / p
        r = r / p
        IF (notlas) goto 225
!     .......... row modification ..........
        DO 200 j = k, n
          p = h (k, j) + q * h (k + 1, j)
          h (k, j) = h (k, j) - p * x
          h (k + 1, j) = h (k + 1, j) - p * y
  200   END DO
!
        j = min0 (en, k + 3)
!     .......... column modification ..........
        DO 210 i = 1, j
          p = x * h (i, k) + y * h (i, k + 1)
          h (i, k) = h (i, k) - p
          h (i, k + 1) = h (i, k + 1) - p * q
  210   END DO
!     .......... accumulate transformations ..........
        DO 220 i = low, igh
          p = x * z (i, k) + y * z (i, k + 1)
          z (i, k) = z (i, k) - p
          z (i, k + 1) = z (i, k + 1) - p * q
  220   END DO
        GOTO 255
  225   CONTINUE
!     .......... row modification ..........
        DO 230 j = k, n
          p = h (k, j) + q * h (k + 1, j) + r * h (k + 2, j)
          h (k, j) = h (k, j) - p * x
          h (k + 1, j) = h (k + 1, j) - p * y
          h (k + 2, j) = h (k + 2, j) - p * zz
  230   END DO
!
        j = min0 (en, k + 3)
!     .......... column modification ..........
        DO 240 i = 1, j
          p = x * h (i, k) + y * h (i, k + 1) + zz * h (i, k + 2)
          h (i, k) = h (i, k) - p
          h (i, k + 1) = h (i, k + 1) - p * q
          h (i, k + 2) = h (i, k + 2) - p * r
  240   END DO
!     .......... accumulate transformations ..........
        DO 250 i = low, igh
          p = x * z (i, k) + y * z (i, k + 1) + zz * z (i, k + 2)
          z (i, k) = z (i, k) - p
          z (i, k + 1) = z (i, k + 1) - p * q
          z (i, k + 2) = z (i, k + 2) - p * r
  250   END DO
  255   CONTINUE
!
  260 END DO
!
      GOTO 70
!     .......... one root found ..........
  270 h (en, en) = x + t
      wr (en) = h (en, en)
      wi (en) = 0.0d0
      en = na
      GOTO 60
!     .......... two roots found ..........
  280 p = (y - x) / 2.0d0
      q = p * p + w
      zz = dsqrt (dabs (q) )
      h (en, en) = x + t
      x = h (en, en)
      h (na, na) = y + t
      IF (q.lt.0.0d0) goto 320
!     .......... real pair ..........
      zz = p + dsign (zz, p)
      wr (na) = x + zz
      wr (en) = wr (na)
      IF (zz.ne.0.0d0) wr (en) = x - w / zz
      wi (na) = 0.0d0
      wi (en) = 0.0d0
      x = h (en, na)
      s = dabs (x) + dabs (zz)
      p = x / s
      q = zz / s
      r = dsqrt (p * p + q * q)
      p = p / r
      q = q / r
!     .......... row modification ..........
      DO 290 j = na, n
        zz = h (na, j)
        h (na, j) = q * zz + p * h (en, j)
        h (en, j) = q * h (en, j) - p * zz
  290 END DO
!     .......... column modification ..........
      DO 300 i = 1, en
        zz = h (i, na)
        h (i, na) = q * zz + p * h (i, en)
        h (i, en) = q * h (i, en) - p * zz
  300 END DO
!     .......... accumulate transformations ..........
      DO 310 i = low, igh
        zz = z (i, na)
        z (i, na) = q * zz + p * z (i, en)
        z (i, en) = q * z (i, en) - p * zz
  310 END DO
!
      GOTO 330
!     .......... complex pair ..........
  320 wr (na) = x + p
      wr (en) = x + p
      wi (na) = zz
      wi (en) = - zz
  330 en = enm2
      GOTO 60
!     .......... all roots found.  backsubstitute to find
!                vectors of upper triangular form ..........
  340 IF (norm.eq.0.0d0) goto 1001
!     .......... for en=n step -1 until 1 do -- ..........
      DO 800 nn = 1, n
        en = n + 1 - nn
        p = wr (en)
        q = wi (en)
        na = en - 1
        IF (q) 710, 600, 800
!     .......... real vector ..........
  600   m = en
        h (en, en) = 1.0d0
        IF (na.eq.0) goto 800
!     .......... for i=en-1 step -1 until 1 do -- ..........
        DO 700 ii = 1, na
          i = en - ii
          w = h (i, i) - p
          r = 0.0d0
!
          DO 610 j = m, en
  610     r = r + h (i, j) * h (j, en)
!
          IF (wi (i) .ge.0.0d0) goto 630
          zz = w
          s = r
          GOTO 700
  630     m = i
          IF (wi (i) .ne.0.0d0) goto 640
          t = w
          IF (t.ne.0.0d0) goto 635
          tst1 = norm
          t = tst1
  632     t = 0.01d0 * t
          tst2 = norm + t
          IF (tst2.gt.tst1) goto 632
  635     h (i, en) = - r / t
          GOTO 680
!     .......... solve real equations ..........
  640     x = h (i, i + 1)
          y = h (i + 1, i)
          q = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i)
          t = (x * s - zz * r) / q
          h (i, en) = t
          IF (dabs (x) .le.dabs (zz) ) goto 650
          h (i + 1, en) = ( - r - w * t) / x
          GOTO 680
  650     h (i + 1, en) = ( - s - y * t) / zz
!
!     .......... overflow control ..........
  680     t = dabs (h (i, en) )
          IF (t.eq.0.0d0) goto 700
          tst1 = t
          tst2 = tst1 + 1.0d0 / tst1
          IF (tst2.gt.tst1) goto 700
          DO 690 j = i, en
            h (j, en) = h (j, en) / t
  690     END DO
!
  700   END DO
!     .......... end real vector ..........
        GOTO 800
!     .......... complex vector ..........
  710   m = na
!     .......... last vector component chosen imaginary so that
!                eigenvector matrix is triangular ..........
        IF (dabs (h (en, na) ) .le.dabs (h (na, en) ) ) goto 720
        h (na, na) = q / h (en, na)
        h (na, en) = - (h (en, en) - p) / h (en, na)
        GOTO 730
  720   CALL cdiv (0.0d0, - h (na, en), h (na, na) - p, q, h (na, na),  &
        h (na, en) )
  730   h (en, na) = 0.0d0
        h (en, en) = 1.0d0
        enm2 = na - 1
        IF (enm2.eq.0) goto 800
!     .......... for i=en-2 step -1 until 1 do -- ..........
        DO 795 ii = 1, enm2
          i = na - ii
          w = h (i, i) - p
          ra = 0.0d0
          sa = 0.0d0
!
          DO 760 j = m, en
            ra = ra + h (i, j) * h (j, na)
            sa = sa + h (i, j) * h (j, en)
  760     END DO
!
          IF (wi (i) .ge.0.0d0) goto 770
          zz = w
          r = ra
          s = sa
          GOTO 795
  770     m = i
          IF (wi (i) .ne.0.0d0) goto 780
          CALL cdiv ( - ra, - sa, w, q, h (i, na), h (i, en) )
          GOTO 790
!     .......... solve complex equations ..........
  780     x = h (i, i + 1)
          y = h (i + 1, i)
          vr = (wr (i) - p) * (wr (i) - p) + wi (i) * wi (i) - q * q
          vi = (wr (i) - p) * 2.0d0 * q
          IF (vr.ne.0.0d0.or.vi.ne.0.0d0) goto 784
          tst1 = norm * (dabs (w) + dabs (q) + dabs (x) + dabs (y)      &
          + dabs (zz) )
          vr = tst1
  783     vr = 0.01d0 * vr
          tst2 = tst1 + vr
          IF (tst2.gt.tst1) goto 783
  784     CALL cdiv (x * r - zz * ra + q * sa, x * s - zz * sa - q * ra,&
          vr, vi, h (i, na), h (i, en) )
          IF (dabs (x) .le.dabs (zz) + dabs (q) ) goto 785
          h (i + 1, na) = ( - ra - w * h (i, na) + q * h (i, en) )      &
          / x
          h (i + 1, en) = ( - sa - w * h (i, en) - q * h (i, na) )      &
          / x
          GOTO 790
  785     CALL cdiv ( - r - y * h (i, na), - s - y * h (i, en), zz, q,  &
          h (i + 1, na), h (i + 1, en) )
!
!     .......... overflow control ..........
  790     t = dmax1 (dabs (h (i, na) ), dabs (h (i, en) ) )
          IF (t.eq.0.0d0) goto 795
          tst1 = t
          tst2 = tst1 + 1.0d0 / tst1
          IF (tst2.gt.tst1) goto 795
          DO 792 j = i, en
            h (j, na) = h (j, na) / t
            h (j, en) = h (j, en) / t
  792     END DO
!
  795   END DO
!     .......... end complex vector ..........
  800 END DO
!     .......... end back substitution.
!                vectors of isolated roots ..........
      DO 840 i = 1, n
        IF (i.ge.low.and.i.le.igh) goto 840
!
        DO 820 j = i, n
  820   z (i, j) = h (i, j)
!
  840 END DO
!     .......... multiply by transformation matrix to give
!                vectors of original full matrix.
!                for j=n step -1 until low do -- ..........
      DO 880 jj = low, n
        j = n + low - jj
        m = min0 (j, igh)
!
        DO 880 i = low, igh
          zz = 0.0d0
!
          DO 860 k = low, m
  860     zz = zz + z (i, k) * h (k, j)
!
          z (i, j) = zz
  880 CONTINUE
!
      GOTO 1001
!     .......... set error -- all eigenvalues have not
!                converged after 30*n iterations ..........
 1000 ierr = en
 1001 RETURN
      END SUBROUTINE hqr2


      SUBROUTINE cdiv (AR, AI, BR, BI, CR, CI)
!     BEGIN PROLOGUE  CDIV
!     SUBSIDIARY
!     PURPOSE  Compute the complex quotient of two complex numbers.
!     LIBRARY   SLATEC
!     TYPE      COMPLEX (CDIV-C)
!     AUTHOR  (UNKNOWN)
!     DESCRIPTION
!
!     Complex division, (CR,CI) = (AR,AI)/(BR,BI)
!
!   SEE ALSO  EISDOC
!   ROUTINES CALLED  (NONE)
!   REVISION HISTORY  (YYMMDD)
!   811101  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   END PROLOGUE  CDIV
      DOUBLE PRECISION :: AR,AI,BR,BI,CR,CI
!
      DOUBLE PRECISION :: S,ARS,AIS,BRS,BIS
!   FIRST EXECUTABLE STATEMENT  CDIV
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END SUBROUTINE cdiv
