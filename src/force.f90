module force
  use global_var, only: wp
  use vector, only: Mol, VecR
  implicit none
! Lennard Jones Potential
! V = 4 * epsilon * [ (sigma/r)**12 - (sigma/r)**6 ]
!   = 4 * epsilon * (sigma/r)**6 * [ (sigma/r)**2 - 1 ]
!   = 4 * r**(-6) * [ r**(-2) - 1 ] for epsilon=sigma=1
! F_i = 48 * epsilon * (sigma/r)**6 * (1/r**2) * [ ( sigma/r)** 6 - 0.5 ] * i where i = x,y,z
!     = 48 * r**(-8) * [ r**(-6) - 0.5 ] * i  for epsilon=sigma=1
  private
  public :: Compute_Force
contains
  subroutine Compute_Force(Atom, cutoff, epsilon, sigma, uSum, region, virsum)
    !! 根据原子的位置和LJ势的参数信息计算原子的受力并转换成加速度, pure过程中无法使用OMP指令对
    type(Mol), intent(inout), dimension(:) :: Atom      !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(in) :: region                    !< 一个二维矢量, 代表盒子在x和y方向上的长度
    real(wp), intent(in) :: cutoff                      !< LJ势的截断半径
    real(wp), intent(in) :: epsilon, sigma              !< LJ势的参数
    real(wp), intent(inout) :: uSum                     !< 模型的总势能
    real(wp), intent(inout) :: virsum            !< 维里应力
    real(wp) :: rCut              !< 截断半径的平方, 这样在比较的时候就不用进行耗费性能的开方处理了
    real(wp) :: rr                !< 原子i和原子j之间距离的平方
    real(wp) :: rri               !< 原子之间的距离的倒数平方
    real(wp) :: rri3              !< 原子之间的距离的倒数的立方, 除法的开销是乘法的好几倍, 应使用中间变量尽量避免除法运算
    real(wp) :: fcval             !< 原子之间的受力
    integer ::  i, j              !< 循环参数
    integer ::  nMol              !< 原子个数
    type(VecR) :: dR              !< 记为第i个原子和第j个原子之间在x和y上的距离

    !! 初始化一些局部变量
    rCut = cutoff * cutoff
    uSum = 0.d0
    Atom%ra = 0.d0    !! 将每个原子的加速度都初始化为0
    virsum = 0.d0
    nMol = size(Atom) !! 计算原子数量

    !! 此处采用verlet算法双层循环更新原子的加速度, 默认x和y均为周期性边界条件, 时间复杂度为N^2
    !$OMP parallel do default(firstprivate) shared(Atom) reduction(+:uSum, virsum) &
    !$OMP SCHEDULE(dynamic, 2 )
    do i = 1, nMol-1
      do j = i + 1, nMol
        dR = Atom(i)%r - Atom(j)%r          !! 注意这里是使用了操作符重载, 才可以直接让2个vecR类的变量相加
        call VWrapAll(dR, region)           !! 考虑周期性边界条件计算原子之间的距离
        rr = dR%VLensq()                    !! 在Fortran中, 幂运算的开销也是乘法开销的好几倍, 应尽量避免类似的运算
        if (rr < rCut) then
          !$OMP CRITICAL
          rri = 1.d0 / rr                   !! 注意这里的rr已经是2个原子之间的距离的平方了
          rri3 = rri * rri * rri            !! 此处是为了避免做除法和幂数运算的写法
          !! 此处可以不需要sigma和epsilon, 但是为了让编译器不报错
          !! 就还是加上了这2个参数, 这2个参数为1,基本不会对仿真有影响
          fcval = 48.d0 * rri3 * (rri3 - 0.5d0) * rri * epsilon * sigma
          call Atom(i)%ra%VSAdd(dR, fcval)    !! 计算第i个原子和第j个原子之间的受力, 根据牛二定律
          call Atom(j)%ra%VSAdd(dR, -fcval)
          uSum = uSum + 4.d0 * rri3 * (rri3 - 1.d0) + 1.d0  !! 计算模型的总势能
          virsum = virsum + fcval * rr                      !! 计算模型的维里应力
          !$OMP END CRITICAL
        end if
      end do
    end do
    !$OMP end parallel do

  end subroutine Compute_Force

  pure subroutine VWrapAll(dR, region)
    !! 考虑周期性边界条件, 计算原子之间的距离
    type(VecR), intent(inout) :: dR
    type(VecR), intent(in) :: region

    if (dR%x >= 0.5 * region%x) then
      dR%x = dR%x - region%x
    else if (dR%x < -0.5 * region%x) then
      dR%x = dR%x + region%x
    end if

    if (dR%y >= 0.5 * region%y) then
      dR%y = dR%y - region%y
    else if (dR%y < -0.5 * region%y) then
      dR%y = dR%y + region%y
    end if

  end subroutine VWrapAll
end module force
