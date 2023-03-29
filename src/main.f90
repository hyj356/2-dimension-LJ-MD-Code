!! 本程序采用最简单的LJ势计算2维情况下的MD动力学
!#define
!#define  FLERR __FILE__,__LINE__       !! __FILE__是预定义字符串宏, 返回其在哪个源文件中, __LINE__是预定义字符串宏, 返回其在哪一行
!! openMP不能用于pure或者elemental的过程中
program main
  use global_var, only: wp, si, stdout, filename, THREADS
  use vector, only: VecR, Mol, vecI, Allocate_Vec
  use force, only: compute_force
  use integerate, only: SingleStep
  use initialize, only: SetParams, SetupJob
  use file_io, only: read_file
  use omp_lib, only: omp_set_num_threads, omp_set_dynamic
  implicit none
  type(Mol), allocatable :: Atoms(:)      !< 一个用于存储所有原子的坐标, 速度, 加速度信息的数组
  type(VecR) :: region, vSum              !< region表示盒子的大小, 而vSum记录整个系统在xy方向上的速度分量
  type(VecI) :: InitUcell                 !< 一个分量是整数的向量, 表示x, y方向上扩展多少个单胞
  real(wp) :: deltaT, density, rCut, temperature, timeNow, uSum, velMag, virSum!, vvSum
  real(wp) :: start, ends
  integer(si) :: moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit, nx, ny, i   !< 仿真用到的参数


  call cpu_time(start)
  !! 从外部文件中读取初始化变量的具体数值
  call read_file(deltaT, density, nx, ny, stepAvg, stepEquil, stepLimit, temperature, InitUcell)

  !! 根据读取到的参数初始化原子的一系列信息
  call SetParams(rCut, region, density, initUcell, nMol, velMag, temperature)

  !! 分配对应的内存
  call Allocate_Vec(Atoms, nMol)

  !! 准备开始计算, 初始化对应的参数
  call SetupJob(Atoms, stepCount, InitUcell, region, vSum, velMag)

  moreCycles = 1
  write(stdout, '(*(A, 5x))') 'step', 'Total', 'KE', 'UE', 'Vsum', 'Pressure'

  !do while(moreCycles == 1)   !! 每次循环stepCount都会加一, 直到达到stepLimit为止, 然后退出循环
  call omp_set_num_threads(THREADS)
  !call omp_set_dynamic(.TRUE.) !! 计算受力的时候, 每一层循环的次数都不一样, 所以动态分配负载更合理
  do i = 1, stepLimit
    call SingleStep(Atoms, stepCount, stepAvg, deltaT, timeNow, region, rCut, uSum, vSum, virSum, density)
    !if (stepCount >= stepLimit ) moreCycles = 0
  end do
  call cpu_time(ends)
  write(stdout, '(A)') "All Computation down!"
  write(stdout, '(A, G0, A)') "本次计算一共耗时: ", (ends - start) / THREADS , ' S'


end program main
