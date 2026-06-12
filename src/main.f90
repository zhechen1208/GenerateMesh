!本程序将用微分方程法计算生成3维O-H网格
!程序由原翅膀网格的一系列生成程序整合修改而成，以期能适应扑旋翼翼的网格生成
!引入了PVM进行并行网格生成，在加快网格生成的同时，能将该部分代码便捷的嵌入3维不可压流的计算程序中
!原翅膀网格系列程序作者：lan等；
!整理修改者：王逗
!日期：2012-3-26
program main
	use mpi	
	use  shared_date_module
	use tip_square_o_module
	implicit none
	double precision x(jmax),y(jmax),z(kmax)   !x，y是等展长平板翼翼型坐标，z是沿展向的坐标分布
	double precision d0,temp,q1,ex,ey,ez,error1,error2,error3,errorth,ex1,ey1,ez1
	double precision d1,thickness,cof,alpha0
	! dr1, dr2 from shared_date_module (set below)
	integer iternum,iter,n1,n2,n3,i1,i2,status(MPI_STATUS_SIZE)
	integer i,j,k,ii,numt,cou
	double precision,external :: stret
	! span_length, domain_length now come from shared_date_module via params.inc
	logical if_plate,goon
	real time_a,time_b,time_c
	real*8 buf1(imax*(jmax+1)),buf2(imax*(jmax+1))
	real s, t
	

	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	print*,"Process",myid,"of",numprocs,"is alive"
	!内存分配
	ipart2=myid*(kmax-2)/numprocs+2
	ipart3=(myid+1)*(kmax-2)/numprocs+1
	k1=ipart2
	k2=ipart3
	if(k1.le.2) k1=2
	if(k2.ge.kmax-1) k2=kmax-1
	print*,myid,k1,k2
	if(myid.gt.0) then
		leftpro=myid-1
	else
		leftpro=MPI_PROC_NULL
	endif
	if(myid.lt.numprocs-1) then
		rightpro=myid+1
	else
		rightpro=MPI_PROC_NULL
	endif

	goon=.false.  !是否属续算

	if(myid.eq.0)then
		open(10,file='output/3d_mesh.dat',status='unknown')
		open(13,file='output/3d_mesh_initial.dat',status='unknown')
		open(12,file='output/goon.dat',status='unknown')
		open(20,file='output/test_error.dat',status='unknown')
	endif
	
	if(goon) then
		if(myid.eq.0)then
			read(12,*)iter
			do k=1,kmax
				DO J=1,Jmax
					DO I=1,Imax
						read(12,'(6f12.6)') X0(i,j,k),Y0(i,j,k),z0(i,j,k),p(i,j,k),q(i,j,k),r(i,j,k)
					enddo
				enddo
			enddo
			close(12)
		endif
		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		call MPI_BCAST(x0,size(x0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(y0,size(y0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(z0,size(z0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(p,size(p),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(q,size(q),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(r,size(r),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		dr1=y0(2,(jmax-1)/4+1,1)-y0(1,(jmax-1)/4+1,1)
		dr2=y0(imax,(jmax-1)/4+1,1)-y0(imax-1,(jmax-1)/4+1,1)
    else	
	
		if_plate=.true.  !机翼是否是平板

		!第一步，获得3维网格的边界
		
!		else !如果是不平板，则用类似翅膀网格系列程序的方法，获得3维网格的边界
		iter=1
		ERR=1.d-6
		pai=3.1415926
		theta0=pai/2.0
		iternum=3000

		thickness=0.01
		alpha0=pai/30.0
		! rad, omega, sigma, dr1_init, dr2_init from shared_date_module / params.inc
		dr1=dr1_init
		dr2=dr2_init
		IF(.not. if_plate) then
		! 	if(myid.eq.0) then
		! 		open(14,file='fore_trailing.dat',status='old')
		! 		open(33,file='test.dat',status='unknown')


		! 		q1=stret(iarc,alpha0,pai/2,alpha) 
		! 		!翼面边界点的生成
		! 		do k=1,kmax
		! 			read(14,*) xi(1,k),yi(1,k),yi((jmax-1)/2+1,k)
		! 			xi((jmax-1)/2+1,k)=xi(1,k)
		! 		enddo

		! 		do i=2,(jmax-1)/2
		! 			do j=1,kmax
		! 				xi(i,j)=xi(1,j)+(xi((jmax-1)/2+1,j)-xi(1,j))/dble((jmax-1)/2)*dble(i-1)
		! 				yi(i,j)=yi(1,j)+(yi((jmax-1)/2+1,j)-yi(1,j))/dble((jmax-1)/2)*dble(i-1)
		! 			enddo
		! 		enddo

		! 		do j=1,kmax
		! 			d1 = thickness*0.5 !*(yi((jmax-1)/2+1,j)-yi(1,j))
		! 			d(1,j)=abs(d1)
		! 			d(2,j)=abs(d1)
		! 		enddo
				
		! 		pi=0.0
		! 		qi=0.0
				
		! 		!----- 迭代生成翼表面网格

		! 		do while(iter.le.inum)
		! 			DO K=1,10000
		! 				ex=0.0
		! 				DO n1=1,3
		! 					CALL XVX1(EX)
		! 				enddo
		! 				ey=0.0
		! 				DO n2=1,3
		! 					CALL YVY1(EY)
		! 				enddo
		! 				!---- 调整左右边界,使其适应内部变化
		! 				do i=1,(jmax-1)/2+1
		! 					yi(i,1)=yi(i,2)
		! 					yi(i,kmax)=yi(i,kmax-1)
		! 				enddo
		! 				if(mod(k,100).eq.0) WRITE(*,'(i6,2f12.8)') K,EX,EY
		! 				ERRxyz=MAX(EX,EY)
		! 				IF(ERRxyz.LE.ERR) exit
		! 			enddo

		! 			if(iter.ne.1) then
		! 				call dpqsur(error1,error2)
		! 				write(*,'(i6,2f12.6)')iter,error1,error2
		! 				if(iter.ge.iternum) exit
		! 				call itpsur
		! 			endif
		! 			iter=iter+1

		! 			!----- 对前后缘作调整
		! 			do j=1,kmax
		! 				!---- 前缘
		! 				i1=1
		! 				i2=iarc
		! 				do i=i1,i2
		! 					cof=alpha(i)
		! 					xi(i,j)=xi(i1,j)+(xi(i2,j)-xi(i1,j))*(1.0-cos(cof))
		! 					yi(i,j)=yi(i1,j)+(yi(i2,j)-yi(i1,j))*(1.0-cos(cof))
		! 				enddo
		! 				!----后缘	
		! 				i1=(jmax-1)/2+1
		! 				i2=(jmax-1)/2-iarc+2
		! 				do i=i1,i2,-1
		! 					cof=alpha((jmax-1)/2+2-i)
		! 					xi(i,j)=xi(i1,j)+(xi(i2,j)-xi(i1,j))*(1.0-cos(cof))
		! 					yi(i,j)=yi(i1,j)+(yi(i2,j)-yi(i1,j))*(1.0-cos(cof))
		! 				enddo
		! 			enddo

		! 			!---- 对左右边界作调整,保证数值边界上展向为直线
		! 			do i=1,(jmax-1)/2+1
		! 				do j=1,kleft-1
		! 					yi(i,j)=yi(i,kleft)
		! 				enddo
		! 			enddo
		! 			do i=1,(jmax-1)/2+1
		! 				do j=kright+1,kmax
		! 					yi(i,j)=yi(i,kright)
		! 				enddo
		! 			enddo
		! 		enddo
		! 		!lai---end while
		! 		!将机翼平面投影投回三维坐标
		! 		do i=1,(jmax-1)/2+1
		! 			ii=(jmax-1)/2+2-i
		! 			do k=1,kmax
		! 				x0(1,i,k)=yi(i,k)
		! 				if(ii.le.iarc) then
		! 					y0(1,i,k)=0.5*thickness*sin(alpha(ii))
		! 				elseif(i.le.iarc) then
		! 					y0(1,i,k)=0.5*thickness*sin(alpha(i))
		! 				else
		! 					y0(1,i,k)=0.5*thickness
		! 				endif
		! 				z0(1,i,k)=xi(i,k)

		! 				x0(1,jmax+1-i,k)=x0(1,i,k)
		! 				y0(1,jmax+1-i,k)=-y0(1,i,k)
		! 				z0(1,jmax+1-i,k)=z0(1,i,k)
		! 			enddo
		! 		enddo
		! 	endif
		else
			if(myid.eq.0)then
				open(11,file='output/innerboundary.dat',status='old')
				read(11,*) 
				DO J=1,kmax
					DO I=1,jmax
						read(11,'(3f12.6)') X0(1,i,j),Y0(1,i,j),z0(1,i,j)
					enddo
				enddo
			endif
		endif

		if(myid.eq.0) then
			print *, '开始为每一个 k 截面生成二维 O 型截面网格...'

			do k = 1, kmax
				call GenerateSectionOMesh(k, rad)
				if (k == 1 .or. k == kmax .or. mod(k,5) == 0) then
					write(*,'(a,i5,a,i5,a,f12.6)') &
						'  已完成截面 k = ', k, ' / ', kmax, &
						', z = ', z0(1,1,k)
				endif
			enddo

			print *, '所有 k 截面二维 O 型网格生成完毕。'
			write(13,*) 'variables = x, y, z'
			write(13,*) "zone i=",imax," ,j=",jmax," ,k=",kmax
			do k=1,kmax
				do j=1,jmax
					do i=1,imax
						write(13,"(3f12.6)") x0(i,j,k), y0(i,j,k), z0(i,j,k)
					enddo
				enddo
			enddo
			close(11)
		endif

		call MPI_BARRIER(MPI_COMM_WORLD,ierr)
		call MPI_BCAST(x0,size(x0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(y0,size(y0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_BCAST(z0,size(z0),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		!初始化源项
		iter=1
		do i=1,imax
			do j=1,jmax
				do k=1,kmax
					p(i,j,k)=0.0
					q(i,j,k)=0.0
					r(i,j,k)=0.0
				enddo
			enddo 
		enddo
		if(myid.eq.0) print*,'边界数据处理完毕！'
	endif

	!计算参数设置
	! num_iter, num_inner, omega, sigma from shared_date_module / params.inc
	eps=0.0001
	ERR=1.d-4
	iternum=300
    pai=3.1415926535897932384626433832795
    theta0=pai/2.0


	if(myid.eq.0) print*,'开始网格点的计算...'
	do while(iter.le.num_iter)
		do K=1,num_inner
			do n1=1,1
				CALL XV1(EX)
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_REDUCE(ex,ex1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) !求各分进程ex总和
				if (myid.eq.0)then
					ex1=sqrt(ex1/((imax-1)*(jmax-1)*(kmax-1)))
				endif
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				numt=imax*(jmax+1)
				!从左向右平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=x0(i,j,k2)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
			   
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.0)then
					cou=0
					do i=1,imax
						do j=0,jmax
							x0(i,j,k1-1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			   
				!从右向左平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=x0(i,j,k1)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.numprocs-1)then
					cou=0
					do i=1,imax
						do j=0,jmax
							x0(i,j,k2+1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			enddo

			DO n2=1,1
				CALL YV1(EY)
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_REDUCE(ey,ey1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) !求各分进程ey总和
				if (myid.eq.0)then
					ey1=sqrt(ey1/((imax-1)*(jmax-1)*(kmax-1)))
				endif
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				numt=imax*(jmax+1)
				!从左向右平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=y0(i,j,k2)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
			   
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.0)then
					cou=0
					do i=1,imax
						do j=0,jmax
							y0(i,j,k1-1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			   
				!从右向左平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=y0(i,j,k1)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.numprocs-1)then
					cou=0
					do i=1,imax
						do j=0,jmax
							y0(i,j,k2+1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			enddo

			do n3=1,1
				call zv1(ez)
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_REDUCE(ez,ez1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr) !求各分进程ex总和
				if (myid.eq.0)then
					ez1=sqrt(ez1/((imax-1)*(jmax-1)*(kmax-1)))
				endif
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				numt=imax*(jmax+1)
				!从左向右平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=z0(i,j,k2)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
			   
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,rightpro,1,buf2,numt,MPI_DOUBLE_PRECISION,leftpro,1,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.0)then
					cou=0
					do i=1,imax
						do j=0,jmax
							z0(i,j,k1-1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			   
				!从右向左平移数据
				cou=0
				do i=1,imax
					do j=0,jmax
						buf1(cou+1)=z0(i,j,k1)
						cou=cou+1
					enddo
				enddo
				call MPI_BARRIER(MPI_COMM_WORLD,ierr)
				call MPI_SENDRECV(buf1,numt,MPI_DOUBLE_PRECISION,leftpro,2,buf2,numt,MPI_DOUBLE_PRECISION,rightpro,2,MPI_COMM_WORLD,status,ierr)
				if(myid.ne.numprocs-1)then
					cou=0
					do i=1,imax
						do j=0,jmax
							z0(i,j,k2+1)=buf2(cou+1)
							cou=cou+1
						enddo
					enddo
				endif
			enddo
			
			if(myid.eq.0) then
				ERRxyz=MAX(EX1,EY1,EZ1)
				if(mod(k,100).eq.0) WRITE(*,'(i6,3f12.6)') K,EX1,EY1,EZ1
			endif
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call MPI_BCAST(ERRxyz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			
			if(iter < 100) then
				err = 5e-5
			elseif(iter < 1000) then
				err = 1e-6
			else
				err = 1e-8
			endif
			IF(ERRxyz.LE.ERR) exit
			ii=k
		enddo
		if(iter.ne.1) then
			call dpqr(error1,error2,error3,errorth)
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call mpi_allreduce(error1,temp,1,MPI_DOUBLE_PRECISION,mpi_max,MPI_COMM_WORLD,ierr)
			error1=temp
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call mpi_allreduce(error2,temp,1,MPI_DOUBLE_PRECISION,mpi_max,MPI_COMM_WORLD,ierr)
			error2=temp
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call mpi_allreduce(error3,temp,1,MPI_DOUBLE_PRECISION,mpi_max,MPI_COMM_WORLD,ierr)
			error3=temp
			if(myid.eq.0)write(*,'(i6,3f12.6)')iter,error1,error2,error3
			e123=max(error1,error2,error3)
			if(e123.le.eps.or.iter.ge.iternum) exit
			!为避免解的震荡，当快收敛时，应适时减小sigma值
			if(error3.lt.100*eps.and.error3.gt.10*eps) then
				!coff=10.0
				coff=1.0
			elseif (error3.lt.10*eps) then
				!coff=10.0
				coff=1.0
			else
				coff=1.0
				!coff=10.0
			endif
			call itp
			if(myid.eq.0) write(20,*) iter,errorth
		endif

		if(iter.gt.100.and.iter.lt.200) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call gather(x0)
			call gather(y0)
			call gather(z0)
			call gather(p)
			call gather(q)
			call gather(r)
			!if(myid.eq.0)then
			!	if(mod(iter,10).eq.0) call out(10)
			!	if(mod(iter,10).eq.0) call out1(iter,12)
			!endif
		endif

		if(iter.gt.200.and.iter.lt.500) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call gather(x0)
			call gather(y0)
			call gather(z0)
			call gather(p)
			call gather(q)
			call gather(r)
			!if(myid.eq.0)then
			!	if(mod(iter,30).eq.0) call out(10)
			!	if(mod(iter,30).eq.0) call out1(iter,12)
			!endif
		endif

		if(iter.gt.500.and.iter.lt.1000) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call gather(x0)
			call gather(y0)
			call gather(z0)
			call gather(p)
			call gather(q)
			call gather(r)
			!if(myid.eq.0)then
			!	if(mod(iter,50).eq.0) call out(10)
			!	if(mod(iter,50).eq.0) call out1(iter,12)
			!	endif
		endif

		if(iter.gt.1000) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call gather(x0)
			call gather(y0)
			call gather(z0)
			call gather(p)
			call gather(q)
			call gather(r)
			!if(myid.eq.0)then
			!	if(mod(iter,100).eq.0) call out(10)
			!	if(mod(iter,100).eq.0) call out1(iter,12)
			!endif
		endif

		if(iter.lt.100) then
			call MPI_BARRIER(MPI_COMM_WORLD,ierr)
			call gather(x0)
			call gather(y0)
			call gather(z0)
			call gather(p)
			call gather(q)
			call gather(r)
			!if(myid.eq.0)then
			!	call out(10)
			!	call out1(iter,12)
			!endif
		endif

		iter=iter+1
	enddo
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call gather(x0)
	call gather(y0)
	call gather(z0)
	call gather(p)
	call gather(q)
	call gather(r)
	if(myid.eq.0) then 
		call out(10)
		call out1(iter,12)

	print *, '开始生成左端block...'
	call GenerateTipBlock(kleft, 1, -1, 'output/left_tip_block.dat')
	print *, '开始生成右端block...'
	call GenerateTipBlock(kright, kmax, 1, 'output/right_tip_block.dat')
		print *, '过渡段网格block生成完毕！'
	endif
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_FINALIZE(ierr)
end program main

subroutine out(unit)
	use  shared_date_module
	integer unit

    write(unit,*) 'variables = x, y, z'
	write(unit,*) 'zone,i=',imax,',j=',jmax,',k=',kmax
	do k=1,kmax
		DO J=1,Jmax
			DO I=1,Imax
				WRITE(unit,'(3f12.6)') X0(i,j,k),Y0(i,j,k),z0(i,j,k)
			enddo
		enddo
	enddo
	rewind(unit)
	return
END
      
      
subroutine out1(iter,unit1)
	use  shared_date_module
	integer iter,unit1
	write(unit1,*)iter
	do k=1,kmax
		DO J=1,Jmax
			DO I=1,Imax
				WRITE(unit1,'(6f12.6)') X0(i,j,k),Y0(i,j,k),z0(i,j,k),p(i,j,k),q(i,j,k),r(i,j,k)
			enddo
		enddo
	enddo
	rewind(unit1)
	return
END

!*******************************************************************
!  COMPUTING X VALUE SUBPROGAM --------- SOR METHOD
!*******************************************************************
SUBROUTINE XV1(EX)
	use  shared_date_module
	implicit none

	double precision EX,XOLD,XE,YE,ZE,XW,YW,ZW,XN,YN,ZN,XS,YS,ZS,XF,YF,ZF,XB,YB,ZB,XEW,YEW,ZEW,XNS,YNS,ZNS,XFB,YFB,ZFB
	double precision T1,T2,T3,TEMP1,TEMP2,TEMP3,A1,A2,A3,B1,B2,B3,COF
	integer i,j,k
	
	x0(:,0,:) = x0(:,Jmax-1,:)
    x0(:,jmax,:) = x0(:,1,:)
	EX=0.0
	do k=k1,k2
		DO I=2,imax-1
			DO J=1,jmax-1
				XOLD=X0(I,J,k)
				XE=X0(I+1,J,k)
				yE=y0(I+1,J,k)
				zE=z0(I+1,J,k)
				xw=x0(i-1,j,k)
				yw=y0(i-1,j,k)
				zw=z0(i-1,j,k)
				xn=x0(i,j+1,k)
				yn=y0(i,j+1,k)
				zn=z0(i,j+1,k)
				xs=x0(i,j-1,k)
				ys=y0(i,j-1,k)
				zs=z0(i,j-1,k)
				xf=x0(i,j,k+1)
				yf=y0(i,j,k+1)
				zf=z0(i,j,k+1)
				xb=x0(i,j,k-1)
				yb=y0(i,j,k-1)
				zb=z0(i,j,k-1)
				XEW=0.5*(XE-XW)
				yEW=0.5*(yE-yW)
				zEW=0.5*(zE-zW)
				XNS=0.5*(XN-XS)
				yNS=0.5*(yN-yS)
				zNS=0.5*(zN-zS)
				xfb=0.5*(xf-xb)
				yfb=0.5*(yf-yb)
				zfb=0.5*(zf-zb)
				t1=xns**2+yns**2+zns**2
				t2=xfb**2+yfb**2+zfb**2
				t3=xew**2+yew**2+zew**2
				temp1=xns*xfb+yns*yfb+zns*zfb
				temp2=xfb*xew+yfb*yew+zfb*zew
				temp3=xew*xns+yew*yns+zew*zns
				a1=t1*t2-temp1**2
				a2=t2*t3-temp2**2
				a3=t3*t1-temp3**2
				b1=temp2*temp1-temp3*t2
				b2=temp3*temp2-temp1*t3
				b3=temp1*temp3-temp2*t1
				cof=OMEGA/(2.0*(a1+a2+a3))
				X0(I,J,k)=(1.0-OMEGA)*X0(I,J,k)+cof*(a1*(XE+XW+p(i,j,k)*xew) &
							+a2*(XN+XS+q(I,J,k)*xns)+a3*(xf+xb+r(i,j,k)*xfb) &
							+0.5*b1*(x0(i+1,j+1,k)-x0(i+1,j-1,k)-x0(i-1,j+1,k)+x0(i-1,j-1,k)) &
							+0.5*b2*(x0(i,j+1,k+1)-x0(i,j+1,k-1)-x0(i,j-1,k+1)+x0(i,j-1,k-1)) &
							+0.5*b3*(x0(i+1,j,k+1)-x0(i-1,j,k+1)-x0(i+1,j,k-1)+x0(i-1,j,k-1)))
				EX=EX+(X0(I,J,k)-XOLD)**2
			enddo
		enddo
	enddo
!	EX=SQRT(EX/((imax-1)*(jmax-1)*(k1-k2)))

	RETURN
END
!*******************************************************************
!  COMPUTING Y VALUE SUBPROGAM --------- SOR METHOD
!*******************************************************************
SUBROUTINE YV1(EY)
	use  shared_date_module
	implicit none

	double precision EY,YOLD,XE,YE,ZE,XW,YW,ZW,XN,YN,ZN,XS,YS,ZS,XF,YF,ZF,XB,YB,ZB,XEW,YEW,ZEW,XNS,YNS,ZNS,XFB,YFB,ZFB
	double precision T1,T2,T3,TEMP1,TEMP2,TEMP3,A1,A2,A3,B1,B2,B3,COF
	integer i,j,k

	y0(:,0,:) = y0(:,Jmax-1,:)
    y0(:,jmax,:) =y0(:,1,:)
	EY=0.0
	do k=k1,k2
		DO I=2,imax-1
			DO J=1,jmax-1
				yOLD=y0(I,J,k)
				XE=X0(I+1,J,k)
				yE=y0(I+1,J,k)
				zE=z0(I+1,J,k)
				xw=x0(i-1,j,k)
				yw=y0(i-1,j,k)
				zw=z0(i-1,j,k)
				xn=x0(i,j+1,k)
				yn=y0(i,j+1,k)
				zn=z0(i,j+1,k)
				xs=x0(i,j-1,k)
				ys=y0(i,j-1,k)
				zs=z0(i,j-1,k)
				xf=x0(i,j,k+1)
				yf=y0(i,j,k+1)
				zf=z0(i,j,k+1)
				xb=x0(i,j,k-1)
				yb=y0(i,j,k-1)
				zb=z0(i,j,k-1)
				XEW=0.5*(XE-XW)
				yEW=0.5*(yE-yW)
				zEW=0.5*(zE-zW)
				XNS=0.5*(XN-XS)
				yNS=0.5*(yN-yS)
				zNS=0.5*(zN-zS)
				xfb=0.5*(xf-xb)
				yfb=0.5*(yf-yb)
				zfb=0.5*(zf-zb)
				t1=xns**2+yns**2+zns**2
				t2=xfb**2+yfb**2+zfb**2
				t3=xew**2+yew**2+zew**2
				temp1=xns*xfb+yns*yfb+zns*zfb
				temp2=xfb*xew+yfb*yew+zfb*zew
				temp3=xew*xns+yew*yns+zew*zns
				a1=t1*t2-temp1**2
				a2=t2*t3-temp2**2
				a3=t3*t1-temp3**2
				b1=temp2*temp1-temp3*t2
				b2=temp3*temp2-temp1*t3
				b3=temp1*temp3-temp2*t1
				cof=OMEGA/(2.0*(a1+a2+a3))
				y0(I,J,k)=(1.0-OMEGA)*y0(I,J,k)+cof*(a1*(yE+yW+p(i,j,k)*yew) &
							+a2*(yN+yS+q(I,J,k)*yns)+a3*(yf+yb+r(i,j,k)*yfb) &
							+0.5*b1*(y0(i+1,j+1,k)-y0(i+1,j-1,k)-y0(i-1,j+1,k)+y0(i-1,j-1,k)) &
							+0.5*b2*(y0(i,j+1,k+1)-y0(i,j+1,k-1)-y0(i,j-1,k+1)+y0(i,j-1,k-1)) &
							+0.5*b3*(y0(i+1,j,k+1)-y0(i-1,j,k+1)-y0(i+1,j,k-1)+y0(i-1,j,k-1)))
				EY=EY+(Y0(I,J,k)-YOLD)**2
			enddo
		enddo
	enddo
!	EY=SQRT(EY/((imax-1)*(jmax-1)*(k1-k2)))
	RETURN
END
!*******************************************************************
!  COMPUTING z VALUE SUBPROGAM --------- SOR METHOD
!*******************************************************************
SUBROUTINE zV1(Ez)
	use  shared_date_module
	implicit none

	double precision EZ,ZOLD,XE,YE,ZE,XW,YW,ZW,XN,YN,ZN,XS,YS,ZS,XF,YF,ZF,XB,YB,ZB,XEW,YEW,ZEW,XNS,YNS,ZNS,XFB,YFB,ZFB
	double precision T1,T2,T3,TEMP1,TEMP2,TEMP3,A1,A2,A3,B1,B2,B3,COF
	integer i,j,k

	z0(:,0,:) = z0(:,Jmax-1,:)
    z0(:,jmax,:) =z0(:,1,:)
	Ez=0.0
	do k=k1,k2
		DO I=2,imax-1
			DO J=1,jmax-1
				zOLD=z0(I,J,k)
				XE=X0(I+1,J,k)
				yE=y0(I+1,J,k)
				zE=z0(I+1,J,k)
				xw=x0(i-1,j,k)
				yw=y0(i-1,j,k)
				zw=z0(i-1,j,k)
				xn=x0(i,j+1,k)
				yn=y0(i,j+1,k)
				zn=z0(i,j+1,k)
				xs=x0(i,j-1,k)
				ys=y0(i,j-1,k)
				zs=z0(i,j-1,k)
				xf=x0(i,j,k+1)
				yf=y0(i,j,k+1)
				zf=z0(i,j,k+1)
				xb=x0(i,j,k-1)
				yb=y0(i,j,k-1)
				zb=z0(i,j,k-1)
				XEW=0.5*(XE-XW)
				yEW=0.5*(yE-yW)
				zEW=0.5*(zE-zW)
				XNS=0.5*(XN-XS)
				yNS=0.5*(yN-yS)
				zNS=0.5*(zN-zS)
				xfb=0.5*(xf-xb)
				yfb=0.5*(yf-yb)
				zfb=0.5*(zf-zb)
				t1=xns**2+yns**2+zns**2
				t2=xfb**2+yfb**2+zfb**2
				t3=xew**2+yew**2+zew**2
				temp1=xns*xfb+yns*yfb+zns*zfb
				temp2=xfb*xew+yfb*yew+zfb*zew
				temp3=xew*xns+yew*yns+zew*zns
				a1=t1*t2-temp1**2
				a2=t2*t3-temp2**2
				a3=t3*t1-temp3**2
				b1=temp2*temp1-temp3*t2
				b2=temp3*temp2-temp1*t3
				b3=temp1*temp3-temp2*t1
				cof=OMEGA/(2.0*(a1+a2+a3))
				z0(I,J,k)=(1.0-OMEGA)*z0(I,J,k)+cof*(a1*(zE+zW+p(i,j,k)*zew) &
						+a2*(zN+zS+q(I,J,k)*zns)+a3*(zf+zb+r(i,j,k)*zfb) &
						+0.5*b1*(z0(i+1,j+1,k)-z0(i+1,j-1,k)-z0(i-1,j+1,k)+z0(i-1,j-1,k)) &
						+0.5*b2*(z0(i,j+1,k+1)-z0(i,j+1,k-1)-z0(i,j-1,k+1)+z0(i,j-1,k-1)) &
						+0.5*b3*(z0(i+1,j,k+1)-z0(i-1,j,k+1)-z0(i+1,j,k-1)+z0(i-1,j,k-1)))
				Ez=Ez+(z0(I,J,k)-zOLD)**2
			enddo
		enddo
	enddo
!	Ez=SQRT(Ez/((imax-1)*(jmax-1)*(k1-k2)))
	RETURN
END

!*******************************************************************
!  COMPUTING X VALUE SUBPROGAM --------- SOR METHOD
!*******************************************************************
SUBROUTINE XVX1(EX)
	use  shared_date_module
	implicit none
	double precision EX,XOLD,XE,XW,XN,XS,YE,YW,YN,YS,XEW,XNS,YEW,YNS,a,b,c,cof
	integer i,j,i1,i0
	EX=0.0
	DO I=2,(jmax-1)/2
		DO J=2,kmax-1
			I1=I+1
			I0=I-1
			if(j.eq.kleft.or.j.eq.kright) cycle
			XOLD=Xi(I,J)
			XE=Xi(I1,J)
			XW=Xi(I0,J)
			XN=Xi(I,J+1)
			XS=Xi(I,J-1)
			YE=Yi(I1,J)
			YW=Yi(I0,J)
			YN=Yi(I,J+1)
			YS=Yi(I,J-1)

			XEW=0.5*(XE-XW)
			XNS=0.5*(XN-XS)
			YEW=0.5*(YE-YW)
			YNS=0.5*(YN-YS)
			A=XNS**2+YNS**2
			B=-(XEW*XNS+YEW*YNS)
			C=XEW**2+YEW**2
			IF(A.EQ.0.0 .AND. C.EQ.0.0) cycle
			cof=OMEGA/(2.*(A+C))

			Xi(I,J)=(1.-OMEGA)*Xi(I,J)+cof*(A*(XE+XW+pi(I,J)*xew)+C*(XN+XS+qi(I,J)*xns)+B*0.5*(Xi(I1,J+1)-Xi(I0,J+1)-Xi(I1,J-1)+Xi(I0,J-1)))
			EX=EX+(Xi(I,J)-XOLD)**2
		enddo
	enddo
	EX=SQRT(EX/((jmax-1)/2*(kmax-1)))

	RETURN
END
!*******************************************************************
!  COMPUTING Y VALUE SUBPROGAM --------- SOR METHOD
!*******************************************************************
SUBROUTINE YVY1(EY)
	use  shared_date_module
	implicit none

	double precision EY,YOLD,XE,XW,XN,XS,YE,YW,YN,YS,XEW,XNS,YEW,YNS,a,b,c,cof
	integer i,j,i1,i0
	!----- 左右边界不动
	EY=0.0
	DO I=2,(jmax-1)/2
		DO J=2,kmax-1
			I1=I+1
			I0=I-1
      
			YOLD=Yi(I,J)
			XE=Xi(I1,J)
			XW=Xi(I0,J)
			XN=Xi(I,J+1)
			XS=Xi(I,J-1)
			YE=Yi(I1,J)
			YW=Yi(I0,J)
			YN=Yi(I,J+1)
			YS=Yi(I,J-1)

			XEW=0.5*(XE-XW)
			XNS=0.5*(XN-XS)
			YEW=0.5*(YE-YW)
			YNS=0.5*(YN-YS)
			A=XNS**2+YNS**2
			B=-(XEW*XNS+YEW*YNS)
			C=XEW**2+YEW**2
			IF(A.EQ.0.0.AND.C.EQ.0.0) cycle

			cof=OMEGA/(2.*(A+C))

			yi(I,J)=(1.-OMEGA)*yi(I,J)+cof*(A*(yE+yW+pi(I,J)*yew)+C*(yN+yS+qi(I,J)*yns)+B*0.5*(yi(I1,J+1)-yi(I0,J+1)-yi(I1,J-1)+yi(I0,J-1)))
			EY=EY+(Yi(I,J)-YOLD)**2
		enddo
	enddo
	EY=SQRT(EY/((jmax-1)/2*(kmax-1)))
	RETURN
END

	subroutine XVxx1(EX)
	  !parameter (II=50,JJ=53)  
	  use shared_date_module
	  implicit none
	  integer i,j,i0,i1 
	  double precision sum,t1,temp,ex,ey,XOLD,XE,XW,XN,XS,YE,YW,YN,YS,XEW,XNS,YEW,YNS,A,B,C,cof
	  sum=0.0
	  t1=0.0
	  do i=2,imax
		temp=sqrt((xa(i,2)-xa(i-1,2))**2+(ya(i,2)-ya(i-1,2))**2)
		sum=sum+temp
	  enddo
	  do i=2,imax-1
		temp=sqrt((xa(i,2)-xa(i-1,2))**2+(ya(i,2)-ya(i-1,2))**2)
		t1=t1+temp
		xa(i,1)=t1/sum*(Xa(imax,1)-Xa(1,1))+Xa(1,1)
		xa(i,jmax)=xa(i,1)
	  enddo

	  EX=0.0
	  do I=2,Imax-1
		  do J=2,Jmax-1
			I1=I+1
			I0=I-1

			XOLD=Xa(I,J)
			XE=Xa(I1,J)
			XW=Xa(I0,J)
			XN=Xa(I,J+1)
			XS=Xa(I,J-1)
			YE=Ya(I1,J)
			YW=Ya(I0,J)
			YN=Ya(I,J+1)
			YS=Ya(I,J-1)

			XEW=0.5*(XE-XW)
			XNS=0.5*(XN-XS)
			YEW=0.5*(YE-YW)
			YNS=0.5*(YN-YS)
			A=XNS**2+YNS**2
			B=-(XEW*XNS+YEW*YNS)
			C=XEW**2+YEW**2
		  IF((A==0.0) .AND. (C==0.0)) then
			EX=EX+(Xa(I,J)-XOLD)**2
		  else
			cof=OMEGA/(2.*(A+C))
			Xa(I,J)=(1.-OMEGA)*Xa(I,J)+cof*(A*(XE+XW+pii(I,J)*xew)+C*(XN+XS+qii(I,J)*xns)+B*0.5*(Xa(I1,J+1)-Xa(I0,J+1)-Xa(I1,J-1)+Xa(I0,J-1)))
			EX=EX+(Xa(I,J)-XOLD)**2
		  endif
		enddo
	 enddo
	 EX=SQRT(EX/((Imax-1)*(Jmax-1)))
	 return
	end
	!计算y坐标值的子过程
	SUBROUTINE YVyy1(EY)  
	  use shared_date_module
	  implicit none
	  integer i,j,i0,i1
	  double precision ex,ey,YOLD,XE,XW,XN,XS,YE,YW,YN,YS,XEW,XNS,YEW,YNS,A,B,C,cof
	  do i=2,imax-1
		ya(i,1)=0.0
		ya(i,jmax)=0.0
	  enddo

	  EY=0.0
	  do I=2,Imax-1
		do J=2,Jmax-1
		  I1=I+1
		  I0=I-1

		  YOLD=Ya(I,J)
		  XE=Xa(I1,J)
		  XW=Xa(I0,J)
		  XN=Xa(I,J+1)
		  XS=Xa(I,J-1)
		  YE=Ya(I1,J)
		  YW=Ya(I0,J)
		  YN=Ya(I,J+1)
		  YS=Ya(I,J-1)

		  XEW=0.5*(XE-XW)
		  XNS=0.5*(XN-XS)
		  YEW=0.5*(YE-YW)
		  YNS=0.5*(YN-YS)
		  A=XNS**2+YNS**2
		  B=-(XEW*XNS+YEW*YNS)
		  C=XEW**2+YEW**2
		  IF(A.EQ.0.0.AND.C.EQ.0.0) then
			EY=EY+(Ya(I,J)-YOLD)**2
		  else
			cof=OMEGA/(2.*(A+C))
			ya(I,J)=(1.-OMEGA)*ya(I,J)+cof*(A*(yE+yW+pii(I,J)*yew)+C*(yN+yS+qii(I,J)*yns)+B*0.5*(ya(I1,J+1)-ya(I0,J+1)-ya(I1,J-1)+ya(I0,J-1)))
			EY=EY+(Ya(I,J)-YOLD)**2
		  endif
		enddo
	 enddo
	 EY=SQRT(EY/((Imax-1)*(Jmax-1)))
	 RETURN
	end

subroutine dpqr(error1,error2,error3,errorth)
	use  shared_date_module
	implicit none

	double precision error1,error2,error3,b,XEW,YEW,ZEW,XFB,YFB,ZFB,XNS,YNS,ZNS,sa,sb,sc,THETA1,THETA2,errorth
	double precision DELQ(imax,jmax,kmax),DELR(imax,jmax,kmax),DELP(imax,jmax,kmax)
	integer i,j,k,ii,jj,kk

	error1=0.0
	error2=0.0
	error3=0.0
	errorth=0.0
	do k=k1,k2
		do j=1,jmax-1

			b=sqrt((x0(2,j,k)-x0(1,j,k))**2+(y0(2,j,k)-y0(1,j,k))**2+(z0(2,j,k)-z0(1,j,k))**2)

			XEW=x0(1,j+1,k)-x0(1,j-1,k)
			yEW=y0(1,j+1,k)-y0(1,j-1,k)
			zEW=z0(1,j+1,k)-z0(1,j-1,k)
			xFB=x0(1,j,k+1)-x0(1,j,k-1)
			yFB=y0(1,j,k+1)-y0(1,j,k-1)
			zFB=z0(1,j,k+1)-z0(1,j,k-1)
			xNS=x0(2,j,k)-x0(1,j,k)
			yNS=y0(2,j,k)-y0(1,j,k)
			zNS=z0(2,j,k)-z0(1,j,k)
			sa = sqrt(xew**2+yew**2+zew**2)
			sb = sqrt(xfb**2+yfb**2+zfb**2)
			sc = sqrt(xns**2+yns**2+zns**2)
			theta1 = acos((xew*xns+yew*yns+zew*zns)/(sa*sc))
			theta2 = acos((xfb*xns+yfb*yns+zfb*zns)/(sb*sc))

			delq(1,j,k)=-sigma*tanh((theta0-theta1)/theta0)
			delr(1,j,k)=-coff*sigma*tanh((theta0-theta2)/theta0)
			delp(1,j,k)=sigma*tanh((dr1-b)/dr1)
			p(1,j,k)=p(1,j,k)+delp(1,j,k)
			q(1,j,k)=q(1,j,k)+delq(1,j,k)
			r(1,j,k)=r(1,j,k)+delr(1,j,k)
			if(abs(delp(1,j,k)).ge.error1) error1=abs(delp(1,j,k))
			if(abs(delq(1,j,k)).ge.error2) error2=abs(delq(1,j,k))
			if(abs(delr(1,j,k)).ge.error3) then 
				error3=abs(delr(1,j,k))
				errorth=theta0-theta2
				ii=1
				jj=j
				kk=k
			endif
			enddo
		enddo

			do k=k1,k2
			do j=1,jmax-1

			b=sqrt((x0(imax-1,j,k)-x0(imax,j,k))**2+(y0(imax-1,j,k)-y0(imax,j,k))**2+(z0(imax-1,j,k)-z0(imax,j,k))**2)

			XEW=x0(imax,j+1,k)-x0(imax,j-1,k)
			yEW=y0(imax,j+1,k)-y0(imax,j-1,k)
			zEW=z0(imax,j+1,k)-z0(imax,j-1,k)
			xFB=x0(imax,j,k+1)-x0(imax,j,k-1)
			yFB=y0(imax,j,k+1)-y0(imax,j,k-1)
			zFB=z0(imax,j,k+1)-z0(imax,j,k-1)
			xNS=x0(imax-1,j,k)-x0(imax,j,k)
			yNS=y0(imax-1,j,k)-y0(imax,j,k)
			zNS=z0(imax-1,j,k)-z0(imax,j,k)
			sa = sqrt(xew**2+yew**2+zew**2)
			sb = sqrt(xfb**2+yfb**2+zfb**2)
			sc = sqrt(xns**2+yns**2+zns**2)
			theta1 = acos((xew*xns+yew*yns+zew*zns)/(sa*sc))
			theta2 = acos((xfb*xns+yfb*yns+zfb*zns)/(sb*sc))

			delq(imax,j,k)=-sigma*tanh((theta0-theta1)/theta0)
			delr(imax,j,k)=-coff*sigma*tanh((theta0-theta2)/theta0)
			delp(imax,j,k)=-sigma*tanh((dr2-b)/dr2)
			p(imax,j,k)=p(imax,j,k)+delp(imax,j,k)
			q(imax,j,k)=q(imax,j,k)+delq(imax,j,k)
			r(imax,j,k)=r(imax,j,k)+delr(imax,j,k)
			if(abs(delp(imax,j,k)).ge.error1) error1=abs(delp(imax,j,k))
			if(abs(delq(imax,j,k)).ge.error2) error2=abs(delq(imax,j,k))
			if(abs(delr(imax,j,k)).ge.error3) then 
				error3=abs(delr(imax,j,k))
				errorth=theta0-theta2
				ii=imax
				jj=j
				kk=k
			endif
		enddo
	enddo
!	print*,errorth,ii,jj,kk
	return
end

subroutine dpqsur(error1,error2)
use  shared_date_module
implicit none

double precision error1,error2,s2,XNS,YNS,XEW,YEW,theta
double precision DELP((jmax-1)/2+1,kmax),DELQ((jmax-1)/2+1,kmax)
integer i,j,i1
error1=0.0
error2=0.0

!----计算距离源项  控制1-5层网格距离
!-------- 计算边界1
      i=1
	i1=iarc
      do j=2,kmax-1
      s2=sqrt((xi(i,j)-xi(i1,j))**2+(yi(i,j)-yi(i1,j))**2)
      delp(i,j)=sigma*tanh((d(1,j)-s2)/d(1,j))
      pi(i,j)=pi(i,j)+delp(i,j)
      if(abs(delp(i,j)).ge.error1) error1=abs(delp(i,j))
      enddo
!-------- 计算边界2
      i=(jmax-1)/2+1
	i1=(jmax-1)/2-iarc+2
      do j=2,kmax-1
      s2=sqrt((xi(i,j)-xi(i1,j))**2+(yi(i,j)-yi(i1,j))**2)
      delp(i,j)=-sigma*tanh((d(2,j)-s2)/d(2,j))
      pi(i,j)=pi(i,j)+delp(i,j)
      if(abs(delp(i,j)).ge.error1) error1=abs(delp(i,j))
      enddo
!----计算角度源项
!-------- 计算边界1
      i=1
      do j=2,kmax-1
!      s1=sqrt((x(i,j)-x(i,j+1))**2+(y(i,j)-y(i,j+1))**2)
      s2=sqrt((xi(i,j)-xi(i+1,j))**2+(yi(i,j)-yi(i+1,j))**2)
!      s3=sqrt((x(i+1,j)-x(i,j+1))**2+(y(i+1,j)-y(i,j+1))**2)
!	tem=(s1**2+s2**2-s3**2)/2.0/s1/s2
!      theta=acos(tem)
	XNS = xi(2,j) - xi(1,j)
	YNS = yi(2,j) - yi(1,j)
	XEW = (xi(1,j+1) - xi(1,j-1))
	YEW = (yi(1,j+1) - yi(1,j-1))
	theta = acos( (XNS*XEW+YNS*YEW)/sqrt( (XNS**2+YNS**2)*(XEW**2+YEW**2) ) )
      delq(i,j)=-0.1*sigma*tanh((theta0-theta)/theta0)
      qi(i,j)=qi(i,j)+delq(i,j)
      if(abs(delq(i,j)).ge.error2) error2=abs(delq(i,j))
      enddo
!-------- 计算边界2
i=(jmax-1)/2+1
      do j=2,kmax-1
!      s1=sqrt((x(i,j)-x(i,j+1))**2+(y(i,j)-y(i,j+1))**2)
      s2=sqrt((xi(i,j)-xi(i-1,j))**2+(yi(i,j)-yi(i-1,j))**2)
!      s3=sqrt((x(i-1,j)-x(i,j+1))**2+(y(i-1,j)-y(i,j+1))**2)
!	tem=(s1**2+s2**2-s3**2)/2.0/s1/s2
!      theta=acos(tem)
	XNS = xi((jmax-1)/2,j) - xi((jmax-1)/2+1,j)
	YNS = yi((jmax-1)/2,j) - yi((jmax-1)/2+1,j)
	XEW = (xi((jmax-1)/2+1,j+1) - xi((jmax-1)/2+1,j-1))
	YEW = (yi((jmax-1)/2+1,j+1) - yi((jmax-1)/2+1,j-1))
	theta = acos((XNS*XEW+YNS*YEW)/SQRT((XNS**2+YNS**2)*(XEW**2+YEW**2)))
      delq(i,j)=-0.1*sigma*tanh((theta0-theta)/theta0)
      qi(i,j)=qi(i,j)+delq(i,j)
      if(abs(delq(i,j)).ge.error2) error2=abs(delq(i,j))
      enddo

      return
      end
subroutine dpqbon(error1,error2)  
	use shared_date_module
	implicit none
	integer i,j
	double precision error1,error2,s2,theta
	double precision XNS(2)
	double precision DELP(imax,jmax),delq(imax,jmax)

	error1=0.0
	error2=0.0
	do j=1,jmax-1
		XNS(1) = xa(2,j)-xa(1,j)
		XNS(2) = ya(2,j)-ya(1,j)
		s2 = sqrt(XNS(1)*XNS(1)+XNS(2)*XNS(2))
		XNS = XNS/s2
		theta=acos(XNS(1)*XT(1,j)+XNS(2)*XT(2,j))
		delq(1,j)=-sigma*tanh((theta0-theta)/theta0)
		delp(1,j)=sigma*tanh((dr1-s2)/dr1)
		pii(1,j)=pii(1,j)+delp(1,j)
		qii(1,j)=qii(1,j)+delq(1,j)
		if(abs(delp(1,j)).ge.error1) error1=abs(delp(1,j))
		if(abs(delq(1,j)).ge.error2) error2=abs(delq(1,j))
		ERROR1 = MAX(ERROR1,ABS(DELP(1,j)))
		ERROR2 = MAX(ERROR2,ABS(DELQ(1,j)))
	end do

	do j=1,jmax-1
		XNS(1) = xa(imax-1,j)-xa(imax,j)
		XNS(2) = ya(imax-1,j)-ya(imax,j)
		s2 = sqrt(XNS(1)*XNS(1)+XNS(2)*XNS(2))
		XNS = XNS/s2
		theta=acos(XNS(1)*XXT(1,j)+XNS(2)*XXT(2,j))

		delq(imax,j)=-sigma*tanh((theta0-theta)/theta0)
		delp(imax,j)=-sigma*tanh((dr2-s2)/dr2)
		pii(imax,j)=pii(imax,j)+delp(imax,j)
		qii(imax,j)=qii(imax,j)+delq(imax,j)
		if(abs(delp(imax,j)).ge.error1) error1=abs(delp(imax,j))
		if(abs(delq(imax,j)).ge.error2) error2=abs(delq(imax,j))
		ERROR1 = MAX(ERROR1,ABS(DELP(imax,j)))
		ERROR2 = MAX(ERROR2,ABS(DELQ(imax,j)))
	enddo

	return
end

subroutine itp
	use  shared_date_module
	implicit none

	double precision mm,cof
	integer i,j,k

	mm=8
	do k=k1,k2
		DO I=2,Imax-1
			DO J=1,Jmax
				cof=exp(-dble(i)/dble(mm))
				p(i,j,k)=p(imax,j,k)+(p(1,j,k)-p(imax,j,k))*cof
				q(i,j,k)=q(imax,j,k)+(q(1,j,k)-q(imax,j,k))*cof
				r(i,j,k)=r(imax,j,k)+(r(1,j,k)-r(imax,j,k))*cof
			enddo
		enddo
	enddo
!	print*,'r(imax,1,2)',r(imax,1,2),'r(1,1,2)',r(1,1,2)
	return
end

subroutine itpsur
	use  shared_date_module
	implicit none

	integer i,j,i1,i2
	!--method 1.
	i1=1
	i2=(jmax-1)/2+1
	do j=2,kmax-1
		do i=i1,i2
			pi(i,j)=pi(i1,j)+(pi(i2,j)-pi(i1,j))/dble(i2-i1)*dble(i-i1)
			qi(i,j)=qi(i1,j)+(qi(i2,j)-qi(i1,j))/dble(i2-i1)*dble(i-i1)
		enddo
	enddo 
end

subroutine itpbon  
	use shared_date_module
	implicit none
	integer i,j,mm
	double precision cof

	mm=8.0
	do j=2,jmax-1
		do i=2,imax-1
		cof=exp(-dble(i)/dble(mm))		  
		pii(i,j)=pii(imax,j)+(pii(1,j)-pii(imax,j))*cof
		qii(i,j)=qii(imax,j)+(qii(1,j)-qii(imax,j))*cof
		enddo
	enddo
	return
end

SUBROUTINE INIT_T()
	use shared_date_module
	IMPLICIT NONE
	double precision,dimension(2):: xew
	double precision :: dy
	integer j

	DO j=2,jmax-1
		xew(1) = xa(1,j+1)-xa(1,j-1)
		xew(2) = ya(1,j+1)-ya(1,j-1)
		dy = sqrt(dot_product(xew,xew))
		XT(:,j) = xew/dy
		xew(1) = xa(imax,j+1)-xa(imax,j-1)
		xew(2) = ya(imax,j+1)-ya(imax,j-1)
		dy = sqrt(dot_product(xew,xew))
		XXT(:,j) = xew/dy
	END DO

	xew(1) = xa(1,2)-xa(1,jmax-1)
	xew(2) = ya(1,2)-ya(1,jmax-1)
	dy = sqrt(xew(1)*xew(1)+xew(2)*xew(2))
	XT(1,1) = xew(1)/dy
	XT(2,1) = xew(2)/dy
	XT(1,jmax) = XT(1,1)
	XT(2,jmax) = XT(2,1)

	xew(1) = xa(imax,2)-xa(imax,jmax-1)
	xew(2) = ya(imax,2)-ya(imax,jmax-1)
	dy = sqrt(xew(1)*xew(1)+xew(2)*xew(2))
	XXT(1,1) = xew(1)/dy
	XXT(2,1) = xew(2)/dy
	XXT(1,jmax) = XXT(1,1)
	XXT(2,jmax) = XXT(2,1)

END SUBROUTINE
   

double precision function stret(m,smin,smax,ss)
	integer m,ii
	double precision smin,smax,ss(m),q0,a,b,d,eps,x,temp
	double precision F
	F(x)=x**(m-1)-d*x+d-1
	if(m==2)then
		ss(1)=0.0
		ss(2)=smax
		stret=1.5
		return
    elseif (m==1) then
		ss(1)=smax
		stret=1.0
		return
    elseif(m==3) then
		ss(1)=0.0
		ss(2)=smin
		ss(3)=smax
		if(((smax-smin)/smin)<1) then
			stret=1.5
            return
        else
			stret=(smax-smin)/smin
            return
        endif 
    else
		eps=0.000001
		ss(1)=0.0
		a=1
		b=((smax-smin)/(m-2.0))/smin
		d=smax/smin
		!在[a，b]范围内用二分法求解F(x)=0的解，解即等比数列公比q
		if((F(b)<0).or.(b<1))then
			print*,'没有存在满足条件的公比值，请重新输入！'
			stret=0.0
			return
		else
			do
				temp=F((a+b)/2)
			if(abs(temp)<eps) then
				stret=(a+b)/2
				do ii=2,(m-1)
					ss(ii)=ss(ii-1)+smin*(((a+b)/2)**(ii-2))
				enddo
				ss(m)=smax	
				return
			elseif(temp>0) then
				b=(a+b)/2
			else
				a=(a+b)/2
			endif
			if(abs(b-a)<eps) then
				stret=(a+b)/2
				do ii=2,(m-1)
					ss(ii)=ss(ii-1)+smin*(((a+b)/2)**(ii-2))
				enddo
				ss(m)=smax
				return
			endif
			enddo  		  		 			
		endif
      endif
end function stret

SUBROUTINE gather(tem)
	use mpi
	use shared_date_module
	implicit none

	real*8 tem(imax,0:jmax,kmax)
	integer sendcount,i,j,k,temsend,i2,i3,ii,i4
	integer,allocatable :: recvcounts(:)
	real*8,allocatable :: recvbuf(:)

	allocate (recvcounts(numprocs),displs(numprocs))

	sendcount=(k2-k1+1)*imax*(jmax+1)
	allocate(sendbuf(sendcount),recvbuf(imax*(jmax+1)*kmax))

	displs(1)=0
	do i=1,numprocs
	i2=(i-1)*(kmax-2)/numprocs+2
	i3=i*(kmax-2)/numprocs+1
	if(i2.le.2) i2=2
	if(i3.ge.kmax-1) i3=kmax-1
	recvcounts(i)=(i3-i2+1)*imax*(jmax+1)
	if(i.ne.1) displs(i)=displs(i-1)+recvcounts(i-1)
	enddo

	i2=1
	do k=k1,k2
	do j=0,jmax
	do i=1,imax
	sendbuf(i2)=tem(i,j,k)
	i2=i2+1
	enddo
	enddo
	enddo
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call mpi_allgatherv(sendbuf,sendcount,mpi_real8,recvbuf,recvcounts,displs,mpi_real8,MPI_COMM_WORLD,ierr)
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)

	do i=1,numprocs
	i2=(i-1)*(kmax-2)/numprocs+2
	i3=i*(kmax-2)/numprocs+1
	if(i2.le.2) i2=2
	if(i3.ge.kmax-1) i3=kmax-1
	ii=displs(i)
	do k=i2,i3
	do j=0,jmax
	do i4=1,imax
	tem(i4,j,k)=recvbuf(ii+1)
	ii=ii+1
	enddo
	enddo
	enddo

	enddo
	deallocate (recvcounts,displs,sendbuf,recvbuf)
end

subroutine GenerateSectionOMesh(ksec, rad_in)
    use shared_date_module
    implicit none

    integer, intent(in) :: ksec
    double precision, intent(in) :: rad_in

    integer :: i, j, iter2d, kk, n1, n2
    double precision :: cof, temp, EX, EY, ERRsec
    double precision :: error1, error2

    !------------------------------------------------------------
    ! 1. 当前 k 截面的内边界与外边界
    !------------------------------------------------------------
    do j = 1, jmax
        cof = (dble(j-1) / dble(jmax-1)) * 2.0d0 * pai

        xa(1,j) = x0(1,j,ksec)
        ya(1,j) = y0(1,j,ksec)

        xa(imax,j) = rad_in * cos(pai - cof)
        ya(imax,j) = rad_in * sin(pai - cof)
    enddo

    !------------------------------------------------------------
    ! 2. 初始代数插值
    !------------------------------------------------------------
    do i = 2, imax-1
        temp = dble(i-1) / dble(imax-1)

        do j = 1, jmax
            xa(i,j) = xa(1,j) + (xa(imax,j) - xa(1,j)) * temp
            ya(i,j) = ya(1,j) + (ya(imax,j) - ya(1,j)) * temp
        enddo
    enddo

    !------------------------------------------------------------
    ! 3. 初始化二维边界控制量
    !------------------------------------------------------------
    call init_T()

    do i = 1, imax
        do j = 1, jmax
            pii(i,j) = 0.0d0
            qii(i,j) = 0.0d0
        enddo
    enddo

    !------------------------------------------------------------
    ! 4. 二维椭圆光顺
    ! 可以先用较少迭代测试，确认拓扑正确后再加。
    !------------------------------------------------------------
    do iter2d = 1, 80

        do kk = 1, 500
            do n1 = 1, 3
                call XVxx1(EX)
            enddo

            do n2 = 1, 3
                call YVyy1(EY)
            enddo

            ERRsec = max(EX, EY)

            if (ERRsec <= 1.0d-5) exit
        enddo

        if (iter2d /= 1) then
            call init_T()
            call dpqbon(error1, error2)

            if (max(error1, error2) <= 1.0d-4) exit

            call itpbon
        endif
    enddo

    !------------------------------------------------------------
    ! 5. 写回三维网格当前 k 截面
    !------------------------------------------------------------
    do j = 1, jmax
        do i = 1, imax
            x0(i,j,ksec) = xa(i,j)
            y0(i,j,ksec) = ya(i,j)

            ! 当前截面是 z = const
            z0(i,j,ksec) = z0(1,j,ksec)
        enddo
    enddo

end subroutine GenerateSectionOMesh

subroutine RelaxFirstLayer(beta)
    use shared_date_module
    implicit none

    real*8, intent(in) :: beta

    integer :: j, k, jp, jm
    real*8 :: tx, ty, tl
    real*8 :: nx1, ny1, nx2, ny2
    real*8 :: vx, vy, dot1, dot2
    real*8 :: xloc, yloc

    do k = 1, kmax
        do j = 1, jmax-1

            jp = j + 1
            jm = j - 1

            if (jp > jmax-1) jp = 1
            if (jm < 1)      jm = jmax-1

            tx = x0(1,jp,k) - x0(1,jm,k)
            ty = y0(1,jp,k) - y0(1,jm,k)

            tl = sqrt(tx*tx + ty*ty)
            if (tl <= 1.0d-14) cycle

            tx = tx / tl
            ty = ty / tl

            nx1 =  ty
            ny1 = -tx

            nx2 = -ty
            ny2 =  tx

            vx = x0(2,j,k) - x0(1,j,k)
            vy = y0(2,j,k) - y0(1,j,k)

            dot1 = vx*nx1 + vy*ny1
            dot2 = vx*nx2 + vy*ny2

            if (dot2 > dot1) then
                nx1 = nx2
                ny1 = ny2
            endif

            xloc = x0(1,j,k) + dr1 * nx1
            yloc = y0(1,j,k) + dr1 * ny1

            x0(2,j,k) = (1.0d0-beta)*x0(2,j,k) + beta*xloc
            y0(2,j,k) = (1.0d0-beta)*y0(2,j,k) + beta*yloc
            z0(2,j,k) = z0(1,j,k)

        enddo

        x0(2,jmax,k) = x0(2,1,k)
        y0(2,jmax,k) = y0(2,1,k)
        z0(2,jmax,k) = z0(2,1,k)

        x0(2,0,k) = x0(2,jmax-1,k)
        y0(2,0,k) = y0(2,jmax-1,k)
        z0(2,0,k) = z0(2,jmax-1,k)

    enddo

end subroutine RelaxFirstLayer