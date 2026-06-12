program plate_3d
	implicit none
	integer i,j,tempk,tempk1
	integer jmax,kmax,N,kleft,kright
	real*8 lamda,t1,t2,gama,pi,theta0,q,temp,temp1,temp2
	real*8 guodu,minlong
	real*8 ,allocatable :: x(:,:,:),y(:,:,:),z(:,:,:),z1(:),t(:)
	double precision,external :: stret
	double precision,external :: yangtiao
	integer celue
	
	pi=acos(-1.0d0)
	celue=2 !平板过渡段的生成策略：1.两个相切圆。2.三次样条曲线
	open(1,file='output/innerboundary.dat',status='unknown')
	!平板机翼参数输入
	print*,'本程序完成3维平板网格平板表面的网格点构造，请按指令输入相应的参数！'
	print*,'请输入周向网格点数(注意：实际点数将自动取为离输入点数最近的模4余1的数)：'
	do
		read*,jmax
		if(jmax<=0) then
			print*,'请输入正的点数！'
		elseif(jmax<9) then
			jmax=9
			print*,'调整后的点数为jmax=9'
			exit
		else
			select case (mod(jmax,4))
			case(0)
				jmax=jmax+1
			case(1)
				jmax=jmax
			case(2)
				jmax=jmax-1
			case(3)
				jmax=jmax+2
			case default
				print*,'未知错误！程序将关闭！'
			end select
			print "(1x,'调整后的点数为jmax=',i5)", jmax
			exit
		endif
	enddo
	print*,'请输入展向网格点数(请输入奇数)：'
	do
		read*,kmax
		if(mod(kmax,2).eq.0) then
			print*,'请输入奇数!'
			print*,'请输入展向网格点数(请输入奇数)：'
		else
			exit
		endif
	enddo
	print*,'请输入平板展弦比(需为正数)：'
	do
		read*,lamda
		if(lamda.le.0) then
			print*,'请输入正数!'
			print*,'请输入平板展弦比(需为正数)：'
		else
			exit
		endif
	enddo
	print*,'请输入计算域展向宽度(大于展弦比)：'
	do
		read*,gama
		if(gama.le.lamda) then
			print*,'请输入大于展弦比的值!'
			print*,'请输入计算域展向宽度(大于展弦比)：'
		else
			exit
		endif
	enddo
	print*, '请输入平板的厚度(0.01≤t≤0.06)：'
	do
		read*,t1
		if(t1.lt.0.01.or.t1.gt.0.06) then
			print*,'请输入0.01至0.06的数!'
			print*,'请输入平板的厚度(0.01≤t≤0.06)：'
		else
			exit
		endif
	enddo
	print*, '请输入圆弧上第一二两点与圆心所成的角度θo(规则：θo为几分之π，就输入几，例如，若θo=π/6，就输入6)：'
	do
	read*,N
	if(N<=2) then
		print*,'请输入大于2的数！'
	else
		theta0=pi/dble(N)
		exit
    endif 
	enddo

	allocate(x(1,jmax,kmax),y(1,jmax,kmax),z(1,jmax,kmax),z1(kmax),t(kmax))
	t2=t1 !到端面处收缩的厚度
	!平板左右端点序号
	kleft=(kmax+1)/2-dnint(dble(kmax)/4.0)
	kright=(kmax+1)/2+dnint(dble(kmax)/4.0)
!	kleft=16
!	kright=46
	print*,'kleft=',kleft,'kright=',kright

	!展向坐标分布,及对应的厚度
	if(celue.eq.1)then
		q=stret((kmax+1)/2-kleft+1,0.25*(t1-t2)*theta0,0.5*lamda+(0.25*pi-0.5)*(t1-t2),z1(kleft:(kmax+1)/2))
		do i=kleft,(kmax+1)/2
			if(z1(i).gt.(0.25*pi*(t1-t2))) then
				tempk=i
				exit
			endif
		enddo
		if(abs(z1(tempk)-(0.25*pi*(t1-t2))).ge.abs(z1(tempk-1)-(0.25*pi*(t1-t2)))) tempk=tempk-1
		if(tempk.le.6) tempk=6
		q=stret(tempk-kleft+1,theta0,pi,z1(kleft:tempk))
		do i=kleft,tempk
			if(z1(i).ge.0.25*pi) then 
				tempk1=i
				exit
			endif
		enddo
		if(abs(z1(tempk1)-0.5*pi).gt.abs(z1(tempk1-1)-0.5*pi)) tempk1=tempk1-1
		if(tempk1.le.4) tempk1=4
		q=stret(tempk1-kleft+1,theta0,0.5*pi,z1(kleft:tempk1))

		if(0.5*pi-z1(tempk1-1).gt.(0.5*pi/(tempk-tempk1))) then
			q=stret(tempk-tempk1+1,0.95*(0.5*pi/(tempk-tempk1)),0.5*pi,z1(tempk1:tempk))
		else
			q=stret(tempk-tempk1+1,0.5*pi-z1(tempk1-1),0.5*pi,z1(tempk1:tempk))
		endif
		z1(tempk1)=0.5*pi

		temp1=0.25*(t1-t2)*(0.5*pi-z1(tempk-1))
		temp2=0.25*(t1-t2)*z1(kleft+1)
		do i=kleft,tempk
			if(i.le.tempk1)then
				t(i)=0.5*t2+0.25*(t1-t2)*(1.0-cos(z1(i)))
				z1(i)=0.25*(t1-t2)*sin(z1(i))		
			else
				t(i)=0.5*t2+0.25*(t1-t2)*(1.0+sin(z1(i)))
				z1(i)=0.25*(t1-t2)*(1.0-cos(z1(i)))+0.25*(t1-t2)
			endif
		enddo
		temp=z1(tempk)
		q=stret((kmax+1)/2-tempk+1,temp1,0.5*lamda-0.5*(t1-t2),z1(tempk:(kmax+1)/2))
		do i=tempk,(kmax+1)/2
			z1(i)=temp+z1(i)
			t(i)=t1*0.5
		enddo

		do i=(kmax+1)/2+1,kright
			z1(i)=2.0*z1((kmax+1)/2)-z1((kmax+1)-i)
			t(i)=t((kmax+1)-i)
		enddo
		temp=z1(kright)
		q=stret(kmax-kright+1,temp2,(gama-lamda)/2.0,z1(kright:kmax))
		do i=kright,kmax
			z1(i)=temp+z1(i)
			t(i)=0.5*t2
		enddo
		temp=z1(kleft)
		q=stret(kleft,temp2,(gama-lamda)/2.0,z1(kleft:1:-1))
		do i=1,kleft
			z1(i)=temp-z1(i)
			t(i)=0.5*t2
		enddo
		z1=-z1


	elseif(celue.eq.2)then
		guodu=min(0.5d0,0.4d0*lamda)
		minlong=0.06*(lamda/5.0)
		q=stret((kmax+1)/2-kleft+1,minlong,0.5*lamda,z1(kleft:(kmax+1)/2))
		tempk = (kmax+1)/2
		do i = kleft+1, (kmax+1)/2
			if (z1(i) .ge. guodu) then
				tempk = i
				exit
			endif
		enddo
		if (tempk > kleft) then
			if (abs(z1(tempk)-guodu) .ge. abs(z1(tempk-1)-guodu)) then
				tempk = tempk - 1
			endif
		endif
		if (tempk <= kleft) tempk = kleft + 1
		q = stret(tempk-kleft+1, minlong, guodu, z1(kleft:tempk))
		do i=kleft,tempk
			t(i)=yangtiao(z1(i),0.0d0,0.5*t2,guodu,0.5*t1)
		enddo
		temp=z1(tempk)
		minlong=q*(z1(tempk)-z1(tempk-1))
		
		q=stret((kmax+1)/2-tempk+1,minlong,0.5*lamda-guodu,z1(tempk:(kmax+1)/2))

		do i=tempk,(kmax+1)/2
			z1(i)=temp+z1(i)
			t(i)=0.5*t1
		enddo
		
		do i=(kmax+1)/2+1,kright
			z1(i)=2.0*z1((kmax+1)/2)-z1((kmax+1)-i)
			t(i)=t((kmax+1)-i)
		enddo
		
		temp=z1(kright)
		minlong=0.07*(lamda/5.0)
		q=stret(kmax-kright+1,minlong,(gama-lamda)/2.0,z1(kright:kmax))
		do i=kright,kmax
			z1(i)=temp+z1(i)
			t(i)=0.5*t2
		enddo
		temp=z1(kleft)
		q=stret(kleft,minlong,(gama-lamda)/2.0,z1(kleft:1:-1))
		do i=1,kleft
			z1(i)=temp-z1(i)
			t(i)=0.5*t2
		enddo
		z1=-z1
	endif
	!每个截面处翼型的生成
	do i=1,kmax
		call profile(t(i),t1,theta0,jmax,x(1,:,i),y(1,:,i))
		z(1,:,i)=z1(i)
	enddo

	!输出表面网格		
	write(1,*) 'zone,i=',1,',j=',jmax,',k=',kmax
	DO J=1,kmax
		DO I=1,jmax
			WRITE(1,'(3f12.6)') X(1,i,j),Y(1,i,j),z(1,i,j)
		enddo
	enddo

end program

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

double precision function yangtiao(x0,x1,y1,x2,y2)
	real*8 x0,x1,x2,y1,y2
	real*8 k,a,b

	k=(y2-y1)/(x2-x1)
	a=-2.0*k/((x1-x2)*(x1-x2))
	b=k*(x1+x2)/((x1-x2)*(x1-x2))
	yangtiao=(x0-x1)*k+y1+(a*x0+b)*(x0-x1)*(x0-x2)
	return
end function yangtiao

subroutine profile(t,t1,theta,jmax,x0,y0)
	implicit none
	integer jmax
	real*8 x0(jmax),y0(jmax),t,t1,theta
	real*8 ,allocatable :: s(:)
	double precision,external :: stret
	integer j0,i,j,j1,tempj
	real*8 pi,minlong,maxlong,q
	
	pi=3.1415926d0
	!第一步：计算四分分之一个平板，先将圆弧段当成平板段，用二分法得到初始公比q
	j0=((jmax-1)/4+1)         !计算四分之一平板上的点数
	minlong=(t1/2.0)*theta
	maxlong=(t1/2.0)*0.5*pi+(1.0-t1)/2
	allocate (s(j0))
	q=stret(j0,minlong,maxlong,s)
	!print*,q
	!第二步：确定圆弧段上与平面上分配的点数
	tempj = j0
	do i = 2, j0
		if (s(i) > ((t1/2.0d0)*0.5d0*pi)) then
			tempj = i
			exit
		endif
	enddo

	if (tempj <= 1) then
		j1 = 2
	else
		if (abs(s(tempj-1) - (t1/2.0d0)*0.5d0*pi) >= &
			abs(s(tempj  ) - (t1/2.0d0)*0.5d0*pi)) then
			j1 = tempj
		else
			j1 = tempj - 1
		endif

		if (j1 < 2) j1 = 2
	endif
	!第三步：分别计算圆弧段上与平面段上的点的分布
	q=stret(j1,theta,0.5*pi,s)
	!print*,q
	minlong=(s(j1)-s(j1-1))*(t1/2)*q
	maxlong=(1.0-t1)/2
	q=stret((j0-j1+1),minlong,maxlong,s(j1))
	!print*,q
	!第四步：确定四分之一平板上的坐标
	do i=1,(j1-1)
		x0(i)=-0.5*cos(s(i))*t1-(1.0-t1)/2.0
!		y0(i)=sin(s(i))*t
		y0(i)=t*sqrt(1-cos(s(i))*cos(s(i)))
	enddo
	do i=j1,j0
		x0(i)=-(1.0-t1)/2+s(i)
		y0(i)=t
	enddo
	!第五步：确定全平板上的坐标点
	do i=(j0+1),((jmax-1)/2+1)
		x0(i)=-x0(j0-(i-j0))
		y0(i)=y0(j0-(i-j0))
	enddo
	j=1
	do i=((jmax-1)/2+2),jmax
		j=j+1
		x0(i)=-x0(j)
		y0(i)=-y0(j)
	enddo
 
end