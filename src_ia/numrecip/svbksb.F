c =====================================================================
c
c     SUBROUTINE svbksb
c
c     PURPOSE:  
c
c              solves A cdot X = B for vector X, where A is specified
c              by the arrays u,w,v as returned by svdcmp. m and n are
c              the logical dimensions o a, and will be equal for square
c              matrices. mp an np are the physical dimensions of a.
c              b(1:m) is the input right-hand side. x(1:n) is
c              the output solution vector. No input quantities are
c              destroyed, so the routine may be called sequentially
c              for different b's.
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 56ff
c
c     SUBROUTINES: 
c
c     VARIABLES:
c 
c     IN:      u       : svd result array u(1:mp,1:np)
c              v       : svd result array v(1:np,1:np)
c              w       : svd result array w(1:np)
c              m       : number of values in first dim of u
c              n       : number of values in first dim of v,w and second of u
c              mp      : first declared dimension of u
c              np      : second declared dimension of u and 
c                          first dimension of v,w
c              b       : input array
c  
c     OUT:     x       : output vector
c
c     REVISIONS:
c
c     date     author   remarks
c
c     26. 7.96   ng     first implementation
c 
c =====================================================================
c
      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      integer m,mp,n,np,nmax
      real b(mp),u(mp,np),v(np,np),w(np),x(np)
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      parameter (nmax=500)
      integer i,j,jj
      real s,tmp(nmax)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      do 12 j=1,n
         s=0.
         if(w(j).ne.0.)then
            do 11 i=1,m
               s=s+u(i,j)*b(i)
 11         continue
            s=s/w(j)
         endif
         tmp(j)=s
 12   continue
      do 14 j=1,n
         s=0.
         do 13 jj=1,n
            s=s+v(j,jj)*tmp(jj)
 13      continue
         x(j)=s
 14   continue
c     
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
