module vector
  use global_var, only: wp, stdout, PI
  implicit none

  private
  public :: VecR, Mol, VecI, Allocate_Vec

  type VecR
    real(wp), public :: x  !< 向量的x分量
    real(wp), public :: y  !< 向量的y分量
  contains
    !! 任何一个派生类与任意一个不属于该类的数据类型做四则运算, 都需要设置至少2个函数, 满足相关运算法则
    procedure, private, pass(self) :: vec_add_vec, vec_add_real, real_add_vec
    procedure, private, pass(self) :: vec_sub_vec, vec_sub_real, real_sub_vec
    procedure, private, pass(self) :: vec_pmul_vec, vec_mul_real, real_mult_vec
    procedure, private, pass(self) :: vec_equal_vec, vec_equal_real
    procedure, private, pass(self) ::vecr_Div_vecI, vecr_Div_vecr
    procedure, public, pass(self) :: VSAdd, VLen, fprintf, VMul, VRand, VPro, VLensq

    generic :: operator(+) => vec_add_vec, vec_add_real, real_add_vec
    generic :: operator(-) => vec_sub_vec, vec_sub_real, real_sub_vec
    generic :: operator(*) => vec_pmul_vec, vec_mul_real, real_mult_vec
    generic :: assignment(=) => vec_equal_vec, vec_equal_real
    generic :: VDiv => vecr_Div_vecI, vecr_Div_vecr   !! 写到这里的时候才发现变量类型错误, 然后就紧急写了个重载弥补
  end type VecR

  type Mol
    type(VecR), public :: r      !< 位置向量
    type(VecR), public  :: rv    !< 速度向量
    type(VecR), public  :: ra    !< 加速度向量
  end type Mol

  type VecI
    integer, public :: x
    integer, public :: y
  end type

  !! 逐元函数(elemental), 意味着传入的即使是数组也可以进行相关的运算操作
contains

  subroutine Allocate_Vec(Atom, nMol)
    !! 为派生类型数组分配对应的内存空间
    type(Mol), intent(inout), allocatable :: Atom(:)      !< 传入的原子集合信息, 包括位置, 速度, 加速度
    integer, intent(in) ::  nMol                          !< 原子个数
    integer ::  Alloc_flag
    character(len=256) :: Err_message

    if (nMol <= 0) then
      write(stdout, '(A)') "在分配内存时出现了负数的元素个数!"
      stop
    end if

    if (allocated(Atom)) THEN
      deallocate(Atom, stat=Alloc_flag, ERRMSG=Err_message)
      if (Alloc_flag /= 0) then
        write(stdout, '(A)') '在释放Atom数组内存时出错, 错误信息如下: '
        write(stdout, '(A)') Err_message
        stop
      end if
    END IF

    allocate(Atom(nMol), stat=Alloc_flag, ERRMSG=Err_message)
    if (Alloc_flag /= 0) then
      write(stdout, '(A)') '在分配Atom数组内存时出错, 错误信息如下: '
      write(stdout, '(A)') Err_message
      stop
    end if

  end subroutine Allocate_Vec

  pure elemental function vec_add_vec(self, vec2) result(res)
    !! 矢量加法, 两个矢量相加
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    type(VecR) :: res !< 返回结果同样为矢量

    res%x = self%x + vec2%x
    res%y = self%y + vec2%y

  end function vec_add_vec

  pure elemental function vec_add_real(self, var) result(res)
    !! 矢量加法, 将一个矢量与一个实数相加
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x + var
    res%y = self%y + var

  end function vec_add_real

  pure elemental function real_add_vec(var,self) result(res)
    !! 矢量加法, 将一个矢量与一个实数相加
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x + var
    res%y = self%y + var

  end function real_add_vec

  pure elemental function vec_sub_vec(self, vec2) result(res)
    !! 矢量减法, 两个矢量相减
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    type(VecR) :: res !< 返回结果同样为矢量

    res%x = self%x - vec2%x
    res%y = self%y - vec2%y

  end function vec_sub_vec

  pure elemental function vec_sub_real(self, var) result(res)
    !! 矢量减法, 将一个矢�������������与一个实��相减
    class(VecR), intent(in) :: self  !< 输入的1个矢量类
    real(wp), intent(in) :: var     !< 要���������������上�����实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x - var
    res%y = self%y - var

  end  function vec_sub_real

  pure elemental function real_sub_vec(var, self) result(res)
    !! 矢量减法, ��一个矢量与一个实数相减
    class(VecR), intent(in) :: self  !< 输入��1个矢量类
    real(wp), intent(in) :: var     !< 要���去的实数
    type(VecR) :: res               !< 返回结果同样为矢量

    res%x = self%x - var
    res%y = self%y - var

  end function real_sub_vec

  pure elemental function vec_pmul_vec(self, vec2) result(res)
    !! 矢量的点乘运算法则
    class(VecR), intent(in) :: self, vec2  !< 输入的2个矢量类
    real(wp) :: res !< 返回双精度的结���

    res = self%x * vec2%x + self%y * vec2%y

  end function vec_pmul_vec

  pure elemental function vec_mul_real(self, var) result(res)
    !! 矢量乘����一个实数
    class(VecR), intent(in) :: self  !< 输入的2��矢量类
    real(wp), intent(in) :: var     !< 要加�����������的实数
    type(VecR) :: res !< ��回1个缩放var���的向量

    res%x = self%x * var
    res%y = self%y * var

  end function vec_mul_real

  pure elemental function real_mult_vec(var, self) result(res)
    !! 矢量乘以一个实数
    class(VecR), intent(in) :: self  !< 输入的2个矢量类
    real(wp), intent(in) :: var     !< 要加上的实数
    type(VecR) :: res !< 返回一个缩放之后的矢量

    res%x = self%x * var
    res%y = self%y * var

  end function real_mult_vec

  pure elemental subroutine vec_equal_vec(self, vec)
    !! 矢量等于一个矢量
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    type(VecR), intent(in) :: vec         !< 另一个矢量

    self%x = vec%x
    self%y = vec%y

  end subroutine vec_equal_vec

  pure elemental subroutine vec_equal_real(self, var)
    !! 让矢量等于一个实数
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    real(wp), intent(in) :: var           !< 实数

    self%x = var
    self%y = var

  end subroutine vec_equal_real

  pure subroutine VSAdd(self, vec, var)
    !! 将self向量加上一个放大了var倍的vec向量
    class(VecR), intent(inout) :: self    !< 输入输出矢量
    type(VecR), intent(in) :: vec         !< 被放大的矢量
    real(wp), intent(in) :: var           !< 放大倍数

    self%x = self%x + vec%x * var
    self%y = self%y + vec%y * var

  end subroutine VSAdd

  pure elemental function VLen(self) result(res)
    !! 返回矢量的长度
    class(VecR), intent(in) :: self   !< 输入矢量
    real(wp) :: res                   !< 矢量长度

    res = sqrt(self%x * self%x + self%y * self%y)

  end function VLen

  pure elemental function VLensq(self) result(res)
    !! 返回矢量的长度的平方
    class(VecR), intent(in) :: self   !< 输入矢量
    real(wp) :: res                   !< 矢量长度

    res = self%x * self%x + self%y * self%y

  end function VLensq

  pure elemental subroutine VMul(self, vec1, vec2)
    !! 将vec1和vec2对应的x和y相乘并赋值给self
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1, vec2

    self%x = vec1%x * vec2%x
    self%y = vec1%y * vec2%y

  end  subroutine VMul

  pure elemental subroutine vecr_Div_vecr(self, vec1, vec2)
    !! 将vec1和vec2对应的x和y相除并赋值给self
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1, vec2

    self%x = vec1%x / vec2%x
    self%y = vec1%y / vec2%y

  end  subroutine vecr_Div_vecr

  pure elemental subroutine vecr_Div_vecI(self, vec1, vec2)
    class(VecR), intent(inout) :: self
    type(VecR), intent(in) :: vec1
    type(vecI), intent(in) :: vec2

    self%x = vec1%x / vec2%x
    self%y = vec1%y / vec2%y

  end subroutine vecr_Div_vecI

  impure elemental subroutine VRand(self)
    !! 将vecr类变量中的分量随机赋予一个值
    class(VecR), intent(inout) :: self
    real(wp) :: seed1, seed2

    call random_seed()
    call random_number(seed1)
    seed2 = 2.d0 * PI * seed1
    self%x =  cos(seed2)
    self%y =  sin(seed2)

  end subroutine VRand

  pure elemental function VPro(self) result(res)
    !! 返回向量坐标相乘的结果
    class(VecR), intent(in) :: self
    real(wp) :: res

    res = self%x * self%y
  end function VPro

  subroutine fprintf(self, fileid)
    !! 将向量以格式化的形式输出
    class(VecR), intent(in) :: self
    integer, intent(in), optional :: fileid     !< 可选的通道id, 如果传入的话就可以输出到文件中

    if (.not. present(fileid) ) then
      write(stdout, '(A, F0.7, A, F0.7, A)') '( ', self%x, ', ', self%y, ' )'
    else
      write(fileid, '( I0, 5X, F0.7, 5x, F0.7)')  self%x, self%y
    end if

  end subroutine fprintf

end module vector
