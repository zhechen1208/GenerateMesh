module tip_square_o_module
    use shared_date_module
    implicit none
contains
subroutine GenerateTipBlock(kStart, kEnd, kStep, filename)
    implicit none

    integer, intent(in) :: kStart, kEnd, kStep
    character(len=*), intent(in) :: filename

    integer :: i, j, kk, kPhys
    integer :: nv1, nv2, nk

    real*8 :: xBond(jmax), yBond(jmax)

    real*8, allocatable :: xCap(:,:), yCap(:,:)
    real*8, allocatable :: xBlock(:,:,:), yBlock(:,:,:), zBlock(:,:,:)

    if (jmax < 5) then
        print *, 'GenerateTipBlock error: jmax is too small.'
        stop
    endif

    if (kStep == 0) then
        print *, 'GenerateTipBlock error: kStep cannot be zero.'
        stop
    endif

    if ((kEnd-kStart)*kStep < 0) then
        print *, 'GenerateTipBlock error: inconsistent kStart, kEnd, kStep.'
        print *, 'kStart=', kStart, ' kEnd=', kEnd, ' kStep=', kStep
        stop
    endif

    !------------------------------------------------------------
    ! k 方向层数
    ! 左端: kStart=kleft,  kEnd=1,    kStep=-1
    ! 右端: kStart=kright, kEnd=kmax, kStep= 1
    !------------------------------------------------------------
    nk = abs(kEnd - kStart) + 1
    
    !------------------------------------------------------------
    ! 二维端面网格点数
    ! nv1: 水平方向，取上下边界配对点数
    ! nv2: 厚度方向层数
    !------------------------------------------------------------
    do j = 1, jmax
        yBond(j) = y0(1,j,1)
    enddo
    call DetermineTipCapSize(jmax, yBond, nv1, nv2)
    print *, 'Tip cap size: nv1=', nv1, ' nv2=', nv2
    allocate(xBlock(nv1,nv2,nk))
    allocate(yBlock(nv1,nv2,nk))
    allocate(zBlock(nv1,nv2,nk))

    !------------------------------------------------------------
    ! 关键修改：
    ! 每一个 k 截面都重新读取轮廓，并重新生成二维 cap。
    ! 不再把 kStart 的 cap 直接 extrude。
    !------------------------------------------------------------
    do kk = 1, nk

        kPhys = kStart + (kk-1)*kStep

        if (kPhys < 1 .or. kPhys > kmax) then
            print *, 'GenerateTipBlock error: kPhys out of range.'
            print *, 'kk=', kk, ' kPhys=', kPhys
            stop
        endif

        do j = 1, jmax
            xBond(j) = x0(1,j,kPhys)
            yBond(j) = y0(1,j,kPhys)
        enddo

        call BuildTipCap(jmax, nv1, nv2, xBond, yBond, xCap, yCap)
        do j = 1, nv2
            do i = 1, nv1
                xBlock(i,j,kk) = xCap(i,j)
                yBlock(i,j,kk) = yCap(i,j)
                zBlock(i,j,kk) = z0(1,1,kPhys)
            enddo
        enddo

        if (allocated(xCap)) deallocate(xCap)
        if (allocated(yCap)) deallocate(yCap)

        if (mod(kk,5) == 0 .or. kk == 1 .or. kk == nk) then
            write(*,'(a,i5,a,i5)') '  generated tip cap layer kk=', kk, ' kPhys=', kPhys
        endif

    enddo

    call WriteTip3DTecplot(filename, xBlock, yBlock, zBlock)

    if (allocated(xBlock)) deallocate(xBlock)
    if (allocated(yBlock)) deallocate(yBlock)
    if (allocated(zBlock)) deallocate(zBlock)

end subroutine GenerateTipBlock

subroutine WriteTip3DTecplot(filename, xBlock, yBlock, zBlock)
    implicit none

    character(len=*), intent(in) :: filename
    real*8, intent(in) :: xBlock(:,:,:), yBlock(:,:,:), zBlock(:,:,:)

    integer :: unit
    integer :: i, j, k
    integer :: nv1, nv2, nk

    nv1 = size(xBlock, 1)
    nv2 = size(xBlock, 2)
    nk  = size(xBlock, 3)

    unit = 95
    open(unit, file=filename, status='unknown')

    write(unit,*) 'variables = x, y, z'
    write(unit,*) 'zone t="tip_block", i=', nv1, ', j=', nv2, ', k=', nk, ', f=point'

    do k = 1, nk
        do j = 1, nv2
            do i = 1, nv1
                write(unit,'(3es20.10)') xBlock(i,j,k), yBlock(i,j,k), zBlock(i,j,k)
            enddo
        enddo
    enddo
    close(unit)
end subroutine WriteTip3DTecplot

subroutine BuildTipCap(jmax_in, nv1, nv2, xBond, yBond, xCap, yCap)
    implicit none

    integer, intent(in) :: jmax_in, nv1, nv2
    real*8, intent(in) :: xBond(jmax_in), yBond(jmax_in)
    real*8, allocatable, intent(out) :: xCap(:,:), yCap(:,:)

    integer :: i, m
    real*8 :: ymin, ymax, yy, eta, s
    real*8 :: xLeft, xRight

    if (nv1 < 2 .or. nv2 < 2) then
        print *, 'BuildTipCap error: nv1 and nv2 must be >= 2.'
        stop
    endif

    allocate(xCap(nv1,nv2), yCap(nv1,nv2))

    ! jmax 通常与 j=1 重合，所以统计范围用 1:jmax-1。
    ymin = minval(yBond(1:jmax_in-1))
    ymax = maxval(yBond(1:jmax_in-1))

    if (abs(ymax-ymin) <= 1.0d-14) then
        print *, 'BuildTipCap error: boundary thickness is too small.'
        stop
    endif

    do m = 1, nv2
        eta = dble(m-1) / dble(nv2-1)

        ! V2 方向分布。
        ! 当前使用 cosine 分布，在上下轮廓附近略加密。
        ! 如果想完全均匀，把下一行改为：yy = ymin + (ymax-ymin)*eta
        yy = ymin + (ymax-ymin) * 0.5d0 * (1.0d0 - cos(acos(-1.0d0)*eta))

        call FindHorizontalSpan(jmax_in, xBond, yBond, yy, xLeft, xRight)

        do i = 1, nv1
            s = dble(i-1) / dble(nv1-1)
            xCap(i,m) = (1.0d0-s)*xLeft + s*xRight
            yCap(i,m) = yy
        enddo
    enddo

end subroutine BuildTipCap


subroutine FindHorizontalSpan(n, x, y, yy, xLeft, xRight)
    implicit none

    integer, intent(in) :: n
    real*8, intent(in) :: x(n), y(n), yy
    real*8, intent(out) :: xLeft, xRight

    integer :: j, cnt
    real*8 :: x1, y1, x2, y2, t, xi
    real*8 :: yminSeg, ymaxSeg, tol
    real*8, allocatable :: xs(:)

    allocate(xs(2*n))
    cnt = 0
    tol = 1.0d-12 * max(1.0d0, maxval(abs(y(1:n))))

    do j = 1, n-1
        x1 = x(j)
        y1 = y(j)
        x2 = x(j+1)
        y2 = y(j+1)

        if (abs(y2-y1) <= tol) then
            ! 水平边界段：如果 yy 落在该水平边上，加入两个端点。
            if (abs(yy-y1) <= tol) then
                call AddUniqueX(xs, cnt, 2*n, x1)
                call AddUniqueX(xs, cnt, 2*n, x2)
            endif
        else
            yminSeg = min(y1, y2)
            ymaxSeg = max(y1, y2)

            if (yy >= yminSeg-tol .and. yy <= ymaxSeg+tol) then
                t = (yy-y1) / (y2-y1)
                if (t >= -tol .and. t <= 1.0d0+tol) then
                    xi = x1 + t*(x2-x1)
                    call AddUniqueX(xs, cnt, 2*n, xi)
                endif
            endif
        endif
    enddo

    if (cnt >= 2) then
        xLeft  = minval(xs(1:cnt))
        xRight = maxval(xs(1:cnt))
    else
        ! 兜底：如果某一层没有找到两个交点，则用全局最左/最右点。
        ! 正常闭合轮廓不应该走到这里。
        xLeft  = minval(x(1:n-1))
        xRight = maxval(x(1:n-1))
        print *, 'FindHorizontalSpan warning: less than two intersections at y = ', yy
    endif

    deallocate(xs)

end subroutine FindHorizontalSpan


subroutine AddUniqueX(xs, cnt, nmax, xnew)
    implicit none

    integer, intent(inout) :: cnt
    integer, intent(in) :: nmax
    real*8, intent(inout) :: xs(nmax)
    real*8, intent(in) :: xnew

    integer :: i
    real*8 :: tol

    tol = 1.0d-11 * max(1.0d0, abs(xnew))

    do i = 1, cnt
        if (abs(xs(i)-xnew) <= tol) return
    enddo

    if (cnt < nmax) then
        cnt = cnt + 1
        xs(cnt) = xnew
    else
        print *, 'AddUniqueX warning: intersection buffer is full.'
    endif

end subroutine AddUniqueX

subroutine DetermineTipCapSize(jmax_in, yBond, nv1, nv2)
    implicit none

    integer, intent(in) :: jmax_in
    real*8, intent(in) :: yBond(jmax_in)
    integer, intent(out) :: nv1, nv2

    integer :: j, l, cnt, bestCnt
    real*8 :: tol, ymaxAbs
    real*8 :: yref

    ymaxAbs = maxval(abs(yBond(1:jmax_in-1)))
    tol = 1.0d-10 * max(1.0d0, ymaxAbs)

    !------------------------------------------------------------
    ! nv1:
    ! 找 y > 0 中出现次数最多的一组 y 值。
    ! 对薄平板端面，这一组就是上表面平直段。
    !------------------------------------------------------------
    bestCnt = 0
    yref = 0.0d0

    do j = 1, jmax_in-1

        if (yBond(j) > tol) then

            cnt = 0

            do l = 1, jmax_in-1
                if (abs(yBond(l) - yBond(j)) <= tol) then
                    cnt = cnt + 1
                endif
            enddo

            if (cnt > bestCnt) then
                bestCnt = cnt
                yref = yBond(j)
            endif

        endif

    enddo

    nv1 = bestCnt

    if (nv1 < 2) then
        print *, 'DetermineTipCapSize warning: failed to determine nv1.'
        print *, 'Use fallback nv1=(jmax+1)/2.'
        nv1 = (jmax_in + 1) / 2
    endif

    !------------------------------------------------------------
    ! nv2:
    ! 统计边界里有多少个不同的 y 层。
    ! 这比写死 17 更合理。
    !------------------------------------------------------------
    call CountUniqueYLevels(jmax_in, yBond, tol, nv2)

    if (nv2 < 2) then
        print *, 'DetermineTipCapSize warning: failed to determine nv2.'
        print *, 'Use fallback nv2=17.'
        nv2 = 17
    endif

end subroutine DetermineTipCapSize

subroutine CountUniqueYLevels(jmax_in, yBond, tol, nLevel)
    implicit none

    integer, intent(in) :: jmax_in
    real*8, intent(in) :: yBond(jmax_in)
    real*8, intent(in) :: tol
    integer, intent(out) :: nLevel

    integer :: j, l
    logical :: isNew

    nLevel = 0

    do j = 1, jmax_in-1

        isNew = .true.

        do l = 1, j-1
            if (abs(yBond(j) - yBond(l)) <= tol) then
                isNew = .false.
                exit
            endif
        enddo

        if (isNew) nLevel = nLevel + 1

    enddo

end subroutine CountUniqueYLevels
end module tip_square_o_module
