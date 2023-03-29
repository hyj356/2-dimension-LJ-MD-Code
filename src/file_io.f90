module file_io
  use global_var, only: wp, stdout, filename, si, NDIM
  use vector, only: VecI, Mol, VecR
  implicit none
  private
  public :: read_file, write_property, Dump_xyz
contains
  subroutine read_file(deltaT, density, nx, ny, stepAvg, stepEquil, stepLimit, temperature, InitUcell)
    !! 按照fortran的namelist格式从外部文件读取仿真参数
    integer(si), intent(inout) :: nx, ny                            !< 在x, y方向上扩展的单胞数
    integer(si), intent(inout) :: stepAvg, stepEquil, stepLimit     !< 输出数据步长, 平均数据步长, 一共仿真多少步
    real(wp), intent(inout) :: deltaT                               !< 时间步长
    real(wp), intent(inout) :: density                              !< 原子密度
    real(wp), intent(inout) :: temperature                          !< 期望的温度
    type(VecI), intent(inout) :: InitUcell                          !< 存放读取到的nx和ny
    namelist /InputVar/ deltaT, density, nx, ny, stepAvg, stepEquil, stepLimit, temperature   !< 以namelist的方式读取
    logical :: flag         !< 判断文件是否存在的逻辑变量
    integer :: io_flag      !< 判断读取文件是否成功的整数
    integer :: inputfile    !< 输入文件通道的整数id

    inquire(file=filename, exist=flag)  !! 查询文件是否在当前路径
    if (.not. flag) then
      write(stdout, *) "File 'parameter.nml' doesn't exist in current path!"
      stop
    end if

    open(newunit=inputfile, file=filename, action='read')       !! 打开文件
    read(inputfile, nml=InputVar, iostat=io_flag)               !! 按照格式初始化仿真参数, 采用nml读取参数就必须使用iostat接受报错信息
    close(inputfile)                                            !! 关闭文件
    write(stdout, '(A)') "The variable for simulation is:"      !! 输出读取到的变量的内容
    write(stdout,  nml=InputVar)
    InitUcell%x = nx
    InitUcell%y = ny

  end subroutine read_file

  subroutine write_property(stepCount, Atom, vSum, uSum, virsum, density)
    !! 将计算图中相关的物理量输出到标准输出(stdin)中
    type(Mol), intent(in), dimension(:) :: Atom   !< 传入的原子集合信息, 包括位置, 速度, 加速度
    type(VecR), intent(inout) :: vSum             !< 系统在x和y方向上的速度分量
    real(wp), intent(in) :: uSum                  !< 系统的平均势能
    real(wp), intent(in) :: virsum, density       !< virsum是virial应力中的一项, 而density用于计算真正的virial应力
    integer, intent(in) :: stepCount              !< 当前的步数
    REAL(WP) :: nMoli                                 !< 原子数量的倒数, 为了便于计算
    real(wp) :: kinEnergy                             !< 动能
    real(wp) :: totEnergy                             !< 系统的总能量
    real(wp) :: sumV                                  !< 总速度
    real(wp) :: vvsum                                 !< 所有粒子的速度平方和
    real(wp) :: pressure                              !< 真正的virial应力

    nMoli = 1.d0 / size(Atom)         !! 为了尽量避免除法, 需要多次用到的分母要提前准备好
    vvsum = sum(Atom%rv%VLensq())     !! 求取

    vSum%x = sum(Atom%rv%x);  vSum%y = sum(Atom%rv%y)   !! 计算系统的速度分量
    sumV = (vSum%x + vSum%y) * nMoli                    !! 计算系统的总动量
    kinEnergy = 0.5 * sum(Atom%rv%x * Atom%rv%x + Atom%rv%y * Atom%rv%y) * nMoli  !! 计算平均动能
    totEnergy = kinEnergy + uSum*nMoli                                            !! 计算平均总能量
    pressure = density * (vvsum + virsum) /(size(Atom) * NDIM)
    write(stdout, "(I0, 3x, *(f0.7, 3x))") stepCount, totEnergy, kinEnergy, uSum*nMoli, sumV, pressure

  end subroutine write_property

  subroutine Dump_xyz(Atom, prefix, frame)
    !! 将固定帧数的原子坐标输出到文件当中, 以供后处理
    type(Mol), intent(in), dimension(:) :: Atom   !< 传入的原子集合信息, 包括位置, 速度, 加速度
    !type(VecR), intent(in) :: region             !< 盒子的大小
    character(len=*), intent(in) :: prefix        !< 输出文件的名称
    integer, intent(in) :: frame                  !< 当前帧数
    character(len=:), allocatable :: fname        !< 结合prefix和frame的文件名称
    character(len=30) ::   suffix                 !< 当前帧数对应的字符串变量
    integer :: dumpid                             !< 文件对应的输出通道的整数id
    integer :: nMol                               !< 原子的个数
    integer :: i                                  !< 循环变量

    allocate(character(len=0) :: fname)
    nMol = size(Atom)
    write(suffix, '(G0)') frame           !! 变量类型转换
    fname = prefix // trim(suffix)
    open(newunit=dumpid, file=fname, action='write')    !! 打开文件
    write(dumpid, '(I0)') nMol
    write(dumpid, '(A)') 'xyz格式的文件中第二行为评论行, 可以写任意的字符串在这里, 当然也可以空出这行不写'
    do i = 1, nMol
      write(dumpid, '(A, E24.15, 2X, E24.15)') 'Ar', Atom(i)%r%x, Atom(i)%r%y
    end do
    close(dumpid)

  end subroutine Dump_xyz
end module
