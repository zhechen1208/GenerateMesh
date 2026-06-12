module shared_date_module
	integer imax,jmax,kmax
	integer kleft,kright,iarc,ipart2,ipart3,k1,k2
	integer num_iter,num_inner
	double precision span_length,domain_length
	double precision rad,omega,sigma,dr1,dr2
	double precision dr1_init,dr2_init
	include 'params.inc'
	double precision x0(imax,0:jmax,kmax),y0(imax,0:jmax,kmax),z0(imax,0:jmax,kmax) !网格坐标
	double precision x1(imax,jmax),y1(imax,jmax)  !左边界网格坐标
	double precision x2(imax,jmax),y2(imax,jmax)  !右边界网格坐标
	double precision xa(imax,jmax),ya(imax,jmax)  !临时边界网格坐标
	double precision XT(2,jmax),XXT(2,jmax)       !内外边界切线方向
	double precision xi((jmax-1)/2+1,kmax),yi((jmax-1)/2+1,kmax) !xi，yi是非等展长平板翼表面的平面投影坐标
	double precision d(2,kmax),alpha(iarc)
	double precision p(imax,0:jmax,kmax),q(imax,0:jmax,kmax),r(imax,0:jmax,kmax),pi((jmax-1)/2+1,kmax),qi((jmax-1)/2+1,kmax),pii(imax,jmax),qii(imax,jmax) !源项
	double precision eps,ERR,pai,theta0,ERRXYZ,e123
	double precision coff
	!mpi公共变量
	integer ierr,myid,numprocs,rightpro,leftpro
	integer,allocatable :: datecount(:),displs(:)
	real*8,allocatable :: sendbuf(:),buff(:,:)
end module 