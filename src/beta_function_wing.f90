program beta_function_wing
    implicit none
    integer :: i,j,k,jmax,kmax,N,kleft,kright
    integer :: nWing, ios
    real*8 :: AR, r1hat, r2hat, pBeta, qBeta, betaB
    real*8 :: cbar, R, gama, pi, theta0, thickRatio, fullThick, halfThick
    real*8 :: cminRatio, cmin, zroot, ztip, s, rhat, chord
    real*8, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:), z1(:), csec(:), hsec(:)
    double precision, external :: beta_chord_hat
    logical :: from_file
    character(len=80) :: line

    pi = acos(-1.0d0)
    from_file = .false.

    open(1,file='output/innerboundary.dat',status='unknown')

    ! try reading from beta_params.txt first
    open(2,file='beta_params.txt',status='old',iostat=ios)
    if (ios == 0) then
        read(2,*,iostat=ios) jmax
        read(2,*,iostat=ios) kmax
        read(2,*,iostat=ios) AR
        read(2,*,iostat=ios) r1hat
        read(2,*,iostat=ios) gama
        read(2,*,iostat=ios) thickRatio
        read(2,*,iostat=ios) cminRatio
        read(2,*,iostat=ios) N
        close(2)
        from_file = .true.
        print*, '读取 beta_params.txt 成功'
    else
        print*, 'beta_params.txt 不存在或不完整，切换为交互输入'
    endif

    if (.not. from_file) then
        print*,'本程序生成 beta-function wing 的三维内边界网格 innerboundary.dat'
        print*,'输出格式与 3d_plate.f90 一致：zone,i=1,j=jmax,k=kmax'
        print*,'（也可以创建 beta_params.txt 避免交互输入，一行一个值）'

        print*,'请输入周向网格点数 jmax（自动调整为模4余1）：'
        do
            read*, jmax
            if (jmax <= 0) then
                print*,'请输入正的点数！'
            elseif (jmax < 9) then
                jmax = 9
                print*, '调整后的点数为 jmax=9'
                exit
            else
                select case (mod(jmax,4))
                case(0); jmax = jmax + 1
                case(1); jmax = jmax
                case(2); jmax = jmax - 1
                case(3); jmax = jmax + 2
                end select
                print "(1x,'调整后的点数为 jmax=',i6)", jmax
                exit
            endif
        enddo

        print*,'请输入展向网格点数 kmax（请输入奇数）：'
        do
            read*, kmax
            if (mod(kmax,2) == 0) then
                print*,'请输入奇数！'
            else
                exit
            endif
        enddo

        print*,'请输入翼展比 AR=R/cbar，论文中为 2 到 6：'
        do
            read*, AR
            if (AR <= 0.0d0) then
                print*,'请输入正数！'
            else
                exit
            endif
        enddo

        print*,'请输入 beta wing 的径向质心位置 r1hat，论文中为 0.4, 0.5, 0.6：'
        do
            read*, r1hat
            if (r1hat <= 0.0d0 .or. r1hat >= 1.0d0) then
                print*,'请输入 0 到 1 之间的数！'
            else
                exit
            endif
        enddo

        print*,'请输入计算域展向宽度 gama，需大于 AR：'
        do
            read*, gama
            if (gama <= AR) then
                print*,'请输入大于 AR 的值！'
            else
                exit
            endif
        enddo

        print*,'请输入厚度比例 thickness/cbar，论文取 0.05：'
        do
            read*, thickRatio
            if (thickRatio <= 0.0d0) then
                print*,'请输入正数！'
            else
                exit
            endif
        enddo

        print*,'请输入端部最小弦长比例 cmin/cbar，例如 0.08；若输入 <=0 则自动取 max(0.02,2.5*thickness)：'
        read*, cminRatio

        print*, '请输入前后缘圆角第一层角度参数 N，例如 θ0=π/6 则输入 6：'
        do
            read*, N
            if (N <= 2) then
                print*,'请输入大于 2 的数！'
            else
                theta0 = pi/dble(N)
                exit
            endif
        enddo
    else
        ! validate params read from file
        if (jmax < 9) jmax = 9
        select case (mod(jmax,4))
        case(0); jmax = jmax + 1
        case(2); jmax = jmax - 1
        case(3); jmax = jmax + 2
        end select
        if (mod(kmax,2) == 0) kmax = kmax + 1
        if (AR <= 0.0d0) AR = 4.0d0
        if (r1hat <= 0.0d0 .or. r1hat >= 1.0d0) r1hat = 0.5d0
        if (gama <= AR) gama = AR + 6.0d0
        if (thickRatio <= 0.0d0) thickRatio = 0.05d0
        if (N <= 2) N = 6
        theta0 = pi/dble(N)
    endif

    cbar = 1.0d0
    R = AR*cbar
    fullThick = thickRatio*cbar
    halfThick = 0.5d0*fullThick

    if (cminRatio <= 0.0d0) then
        cminRatio = max(0.02d0, 2.5d0*thickRatio)
    endif
    cmin = cminRatio*cbar

    ! Ellington beta wing parameters.
    r2hat = 0.929d0 * r1hat**0.732d0
    pBeta = r1hat * ( r1hat*(1.0d0-r1hat)/(r2hat*r2hat-r1hat*r1hat) - 1.0d0 )
    qBeta = pBeta * (1.0d0-r1hat)/r1hat
    betaB = exp(log_gamma(pBeta) + log_gamma(qBeta) - log_gamma(pBeta+qBeta))

    print*, 'r2hat = ', r2hat
    print*, 'pBeta = ', pBeta, ' qBeta = ', qBeta
    print*, 'cmin/cbar = ', cminRatio

    allocate(x(1,jmax,kmax), y(1,jmax,kmax), z(1,jmax,kmax))
    allocate(z1(kmax), csec(kmax), hsec(kmax))

    ! 与原 3d_plate.f90 保持一致：机翼主体约占中间一半展向网格。
    kleft  = (kmax+1)/2 - nint(dble(kmax)/4.0d0)
    kright = (kmax+1)/2 + nint(dble(kmax)/4.0d0)
    nWing = kright - kleft + 1

    print*, 'kleft=', kleft, 'kright=', kright

    zroot = -0.5d0*R
    ztip  =  0.5d0*R

    ! 左侧背景段：保持很小弦长，方便原 main.f90 的单 block 生成。
    if (kleft > 1) then
        do k = 1, kleft-1
            s = dble(k-1)/dble(kleft-1)
            z1(k) = -0.5d0*gama + s*(zroot + 0.5d0*gama)
            csec(k) = cmin
            hsec(k) = min(halfThick, 0.45d0*csec(k))
        enddo
    endif

    ! beta-function wing 主体：rhat 从 root 到 tip。
    do k = kleft, kright
        s = dble(k-kleft)/dble(kright-kleft)
        ! 余弦分布，在 root/tip 附近略加密。
        rhat = 0.5d0*(1.0d0 - cos(pi*s))
        z1(k) = zroot + R*rhat
        chord = cbar * beta_chord_hat(rhat, pBeta, qBeta, betaB)
        csec(k) = max(chord, cmin)
        hsec(k) = min(halfThick, 0.45d0*csec(k))
    enddo

    ! 右侧背景段。
    if (kright < kmax) then
        do k = kright+1, kmax
            s = dble(k-kright)/dble(kmax-kright)
            z1(k) = ztip + s*(0.5d0*gama - ztip)
            csec(k) = cmin
            hsec(k) = min(halfThick, 0.45d0*csec(k))
        enddo
    endif

    ! 生成每个 spanwise 截面的闭合薄板轮廓。
    do k = 1, kmax
        call profile_beta(csec(k), hsec(k), theta0, jmax, x(1,:,k), y(1,:,k))
        z(1,:,k) = z1(k)
    enddo

    write(1,*) 'zone,i=',1,',j=',jmax,',k=',kmax
    do k = 1, kmax
        do j = 1, jmax
            write(1,'(3f12.6)') x(1,j,k), y(1,j,k), z(1,j,k)
        enddo
    enddo

    close(1)
    print*, 'innerboundary.dat 已生成。'

end program beta_function_wing


double precision function beta_chord_hat(rhat, p, q, betaB)
    implicit none
    real*8, intent(in) :: rhat, p, q, betaB
    real*8 :: rr, epsr

    epsr = 1.0d-8
    rr = max(epsr, min(1.0d0-epsr, rhat))
    beta_chord_hat = rr**(p-1.0d0) * (1.0d0-rr)**(q-1.0d0) / betaB
end function beta_chord_hat


double precision function stret(m,smin,smax,ss)
    implicit none
    integer, intent(in) :: m
    integer :: ii
    double precision, intent(in) :: smin,smax
    double precision, intent(out) :: ss(m)
    double precision :: a,b,d,eps,x,temp
    double precision :: F

    F(x)=x**(m-1)-d*x+d-1

    if (m == 1) then
        ss(1) = smax
        stret = 1.0d0
        return
    elseif (m == 2) then
        ss(1) = 0.0d0
        ss(2) = smax
        stret = 1.5d0
        return
    elseif (m == 3) then
        ss(1) = 0.0d0
        ss(2) = smin
        ss(3) = smax
        if (((smax-smin)/smin) < 1.0d0) then
            stret = 1.5d0
        else
            stret = (smax-smin)/smin
        endif
        return
    else
        eps = 1.0d-6
        ss(1) = 0.0d0
        a = 1.0d0
        b = ((smax-smin)/dble(m-2))/smin
        d = smax/smin
        if ((F(b) < 0.0d0) .or. (b < 1.0d0)) then
            print*,'stret: 没有存在满足条件的公比值，改用均匀分布。m=',m
            do ii = 1, m
                ss(ii) = smax*dble(ii-1)/dble(m-1)
            enddo
            stret = 1.0d0
            return
        else
            do
                temp = F((a+b)/2.0d0)
                if (abs(temp) < eps) then
                    stret = (a+b)/2.0d0
                    do ii = 2, m-1
                        ss(ii) = ss(ii-1) + smin*stret**(ii-2)
                    enddo
                    ss(m) = smax
                    return
                elseif (temp > 0.0d0) then
                    b = (a+b)/2.0d0
                else
                    a = (a+b)/2.0d0
                endif
                if (abs(b-a) < eps) then
                    stret = (a+b)/2.0d0
                    do ii = 2, m-1
                        ss(ii) = ss(ii-1) + smin*stret**(ii-2)
                    enddo
                    ss(m) = smax
                    return
                endif
            enddo
        endif
    endif
end function stret


subroutine profile_beta(chord, half_thick, theta, jmax, x0, y0)
    implicit none
    integer, intent(in) :: jmax
    real*8, intent(in) :: chord, half_thick, theta
    real*8, intent(out) :: x0(jmax), y0(jmax)
    real*8, allocatable :: ss(:)
    double precision, external :: stret
    integer :: j0, i, j, j1, tempj
    real*8 :: pi, minlong, maxlong, qrat, flat, h, c

    pi = acos(-1.0d0)
    c = chord
    h = min(half_thick, 0.45d0*c)
    flat = max(c - 2.0d0*h, 1.0d-8)

    j0 = (jmax-1)/4 + 1
    allocate(ss(j0))

    ! 四分之一轮廓：左端圆角鼻尖 -> 上表面中心线。
    minlong = h*theta
    maxlong = h*0.5d0*pi + 0.5d0*flat
    qrat = stret(j0, minlong, maxlong, ss)

    tempj = j0
    do i = 2, j0
        if (ss(i) > h*0.5d0*pi) then
            tempj = i
            exit
        endif
    enddo

    if (tempj <= 1) then
        j1 = 2
    else
        if (abs(ss(tempj-1) - h*0.5d0*pi) >= abs(ss(tempj) - h*0.5d0*pi)) then
            j1 = tempj
        else
            j1 = tempj - 1
        endif
        if (j1 < 2) j1 = 2
    endif

    qrat = stret(j1, theta, 0.5d0*pi, ss)
    minlong = (ss(j1)-ss(j1-1))*h*qrat
    maxlong = 0.5d0*flat
    qrat = stret(j0-j1+1, minlong, maxlong, ss(j1))

    ! left rounded end, local coordinate initially centred at mid-chord.
    do i = 1, j1-1
        x0(i) = -0.5d0*flat - h*cos(ss(i))
        y0(i) =  h*sin(ss(i))
    enddo

    ! upper flat part to mid-chord.
    do i = j1, j0
        x0(i) = -0.5d0*flat + ss(i)
        y0(i) = h
    enddo

    ! mirror to upper right.
    do i = j0+1, (jmax-1)/2 + 1
        x0(i) = -x0(j0-(i-j0))
        y0(i) =  y0(j0-(i-j0))
    enddo

    ! mirror to lower side.
    j = 1
    do i = (jmax-1)/2 + 2, jmax
        j = j + 1
        x0(i) = -x0(j)
        y0(i) = -y0(j)
    enddo

    ! Move origin from mid-chord to quarter-chord line: x=0 is 1/4 chord from leading edge.
    do i = 1, jmax
        x0(i) = x0(i) + 0.25d0*c
    enddo

    deallocate(ss)
end subroutine profile_beta
