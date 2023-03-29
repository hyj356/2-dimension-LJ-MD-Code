module initialize
  use global_var, only: wp
  use vector, only: Mol, VecR, VecI
  implicit none

  private
  public :: SetParams, SetupJob


contains
  subroutine SetupJob(Atom, stepCount, InitUcell, region, vSum, velMag)
    !! 初始化所有仿真参数
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(inout) :: vSum                 !< 动能之和
    real(wp), intent(inout) :: velMag                 !< 速度幅值
    type(VecI), intent(in) :: InitUcell               !< 在x和y方向上的单胞个数
    type(VecR), intent(in) ::  region                 !< 分别对应一个单胞的xy长度和盒子的xy长度
    integer, intent(out) :: stepCount                 !< 将步数清0

    stepCount = 0
    call InitCoord(Atom, InitUcell, region)
    call InitVels(Atom, vSum, velMag)
    call InitAccels(Atom)
  end subroutine SetupJob

  subroutine InitCoord(Atom, InitUcell, region)
    !! 初始化原子的坐标
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(in) ::  region                 !< 分别对应一个单胞的xy长度和盒子的xy长度
    type(VecI), intent(in) :: InitUcell               !< 在x和y方向上的单胞个数
    integer :: n, nx, ny                              !< 一共有n个单胞, x方向上有nx个, y方向上有ny个, n = nx * ny
    type(VecR) :: c, gap, tempc                       !< 计算用到的临时变量

    call gap%VDiv(region, InitUcell)
    n = 1

    !! 双层循环放置原子
    do ny = 0, InitUcell%y - 1
      do nx = 0, InitUcell%x - 1
        c%x = nx + 0.5d0;   c%y = ny + 0.5d0
        tempc = c
        call c%VMul(tempc, gap)
        call c%VSAdd(region, -0.5d0)
        Atom(n)%r = c
        n = n + 1
      end do
    end do

    !! debug使用
    ! open(unit=999, file='coord.txt', action='write')
    ! write(999, "(A, *(G0, 5x))") '# ', 'Atom_id', 'x_coord', 'y_coord'
    ! do nx = 1, 400
    !   write(999, "(*(G0, 5x))") nx, Atom(nx)%r%x, Atom(nx)%r%y
    ! end do
    !! debug通过

  end subroutine InitCoord

  subroutine InitVels(Atom, vSum, velMag)
    !! 初始化所有原子的速度以保证其温度达到目标要求, 同时保证系统质心速度分量为0
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(inout) :: vSum    !< 动能之和
    real(wp), intent(in) :: velMag       !< 速度放缩的系数
    real(wp) :: nMol1                    !< nMol的倒数, 为了加快计算效率, 避免每次循环都算一次除法, 设置的中间变量
    type(VecR) :: nMol1_vSum             !< 同样是为了减少计算量, 设置的一个vecR类的中间变量
    integer :: nMol                      !< 原子个数
    integer :: i                         !< 循环变量

    !! 初始化一些变量
    nMol = size(Atom)
    vSum = 0.d0
    nMol1 = -1.0/nMol

    call Atom%rv%VRand()             !! 随机分配一个0-1之间均匀分布的速度
    Atom%rv = Atom%rv * velMag       !! 进行缩放
    vSum%x = sum(Atom%rv%x)          !! 速度累加
    vSum%y = sum(Atom%rv%y)
    nMol1_vSum = nMol1 * vSum

    do i = 1, nMol
      Atom(i)%rv = Atom(i)%rv + nMol1_vSum    !! 使模型的质心动量为0, 但是当然不可能完全为0 ,只能在有限范围内尽量接近0
    end do

    !! debug使用
    ! open(unit=999, file='vel.txt', action='write')
    ! write(999, "(A, *(G0, 5x))") '# ', 'Atom_id', 'x_vel', 'y_vel'
    ! do i = 1, 400
    !   write(999, "(*(G0, 5x))") i, Atom(i)%rv%x, Atom(i)%rv%y
    ! end do
    ! close(999)
    !! debug通过

  end subroutine InitVels

  pure subroutine InitAccels(Atom)
    !! 这里只是简单的将所有的原子的加速度都初始化为0
    type(Mol), intent(inout), dimension(:) :: Atom    !< 传入的原子集合信息, 包括位置, 速度, 加速度

    Atom%ra = 0.d0

  end subroutine InitAccels

  pure subroutine SetParams(rCut, region, density, initUcell, nMol, velMag, temperature)
    !! 根据从�������部文������读取到的参数计算并初始化参数
    type(VecI), intent(in) :: initUcell
    real(wp), intent(in) :: density, temperature
    integer, intent(out) :: nMol              !! 原子个数
    real(wp), intent(out) :: rCut, velMag     !! 截断半径, 速度
    type(VecR), intent(out) :: region
    integer, parameter :: NDIM = 2            !! 模拟在二���环境下进行

    rCut = 2.0 ** (1./6)                      !! 计算截断半径
    region%x = initUcell%x / sqrt(density)
    region%y = initUcell%y / sqrt(density)
    nMol = initUcell%x * initUcell%y          !! 计算原子个��
    velMag = sqrt(NDIM * (1.d0 - 1.d0 / nMol) * temperature)    !! 根据温度计算公式倒推放缩系数
  end subroutine SetParams

end module initialize
