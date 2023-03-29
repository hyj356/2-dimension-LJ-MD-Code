module integerate
  use global_var, only: wp, stdout
  use vector, only: VecR, Mol
  use force, only: compute_force
  use file_io, only: write_property, Dump_xyz
  implicit none

  private
  public :: SingleStep

contains
  subroutine SingleStep(Atom, stepCount, stepAvg, deltaT, timeNow, region, cutoff, uSum, vSum, virsum, density)
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    real(wp), intent(inout) :: timeNow                !< 当前的时刻
    real(wp), intent(inout) :: uSum                   !< 总势能
    real(wp), intent(inout) :: virsum                 !< 维里应力
    type(VecR), intent(inout) :: vSum                 !< 系统在x和y方向上的速度分量
    integer, intent(inout) :: stepCount               !< 当前的步数
    real(wp), intent(in) :: deltaT                    !< 时间步长
    real(wp), intent(in) :: cutoff                    !< 截断半径
    real(wp), intent(in) :: density                   !< 密度, 用于计算virial应力
    integer, intent(in) :: stepAvg                    !< 平均输出步长
    type(VecR), intent(in) :: region                  !< 当前仿真区域的xy方向的长度
    REAL(WP) :: nMol1                                 !< 原子数量的倒数, 为了便于计算
    real(wp) :: kinEnergy                             !< 动能
    real(wp) :: totEnergy                             !< 系统的总能量
    real(wp) :: sumV                                  !< 总速度


    nMol1 = 1.d0 / size(Atom)             !! 为了尽量避免除法, 需要多次用到的分母要提前准备
    timeNow = stepCount * deltaT          !! 当前的仿真时间
    stepCount = stepCount + 1             !! 仿真步数加一
    if (stepCount == 1 ) then
      uSum = 0.d0; virsum = 0.d0
      call write_property(0, Atom, vSum, uSum, virsum, density)     !! 将相关信息输出到屏幕上面
      call Dump_xyz(Atom, '../xyz/trajectory_p_', 0)
    end if
    call LeapfrogStep(Atom, 1, deltaT)    !! 更新位置
    call ApplyBoundaryCond(Atom, region)  !! 利用周期性边界条件检查原子是否飞出盒子外面

    call Compute_Force(Atom, cutoff, 1.d0, 1.d0, uSum, region, virsum)    !! 计算受力, 更新原子的加速度

    call LeapfrogStep(Atom, 2, deltaT)    !! 更新速度
    if (mod(stepCount, stepAvg) == 0 ) then
      call write_property(stepCount, Atom, vSum, uSum, virsum, density)     !! 将相关信息输出到屏幕上面
      call Dump_xyz(Atom, '../xyz/trajectory_p_', stepCount)
    end if


  end subroutine SingleStep

  pure subroutine LeapfrogStep (Atom, part, deltaT)
    !! 采用蛙跳积分法更新原子的位置和速度
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    real(wp), intent(in) :: deltaT                    !< 时间步长
    integer, intent(in) :: part                       !< 判断变量, 这个变量的值决定更新速度还是加速度
    integer :: i                                      !< 循环变量
    integer :: nMol                                   !< 原子数量

    nMol = size(Atom)

    if (part == 1) then
      !! 更新t + 0.5*deltaT时刻的速度和t + deltaT时刻的位置
      do i = 1, nMol
        call Atom(i)%rv%VSAdd(Atom(i)%ra, 0.5d0 * deltaT)   !! 使用0.5*算法避免进行除法操作
        call Atom(i)%r%VSAdd(Atom(i)%rv, deltaT)
      end do
    else
      !! 更新t + deltaT时刻的速度
      do i = 1, nMol
        call Atom(i)%rv%VSAdd(Atom(i)%ra, 0.5d0 * deltaT)
      end do
    end if

  end subroutine LeapfrogStep

  pure subroutine ApplyBoundaryCond(Atom, region)
    !! 将跑出盒子空间的原子利用周期性边界条件重新wrap回去
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(in) :: region                  !< 当前仿真区域的xy方���的长度
    integer :: nMol, i                                !< nMol为原子数量, i为循������������量

    nMol = size(Atom)

    do i = 1, nMol
      if (Atom(i)%r%x >= 0.5 * region%x) then
        Atom(i)%r%x = Atom(i)%r%x - region%x
      else if (Atom(i)%r%x < -0.5 * region%x) then
        Atom(i)%r%x = Atom(i)%r%x + region%x
      end if

      if (Atom(i)%r%y >= 0.5 * region%y) then
        Atom(i)%r%y = Atom(i)%r%y - region%y
      else if (Atom(i)%r%y < -0.5 * region%y) then
        Atom(i)%r%y = Atom(i)%r%y + region%y
      end if
    end do
  end subroutine ApplyBoundaryCond

end module integerate
