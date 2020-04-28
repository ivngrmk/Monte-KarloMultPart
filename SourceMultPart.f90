program MonteKarlo
    !��, ��� ����� ��� ������ �������.
    use dflib
    type (wxycoord) wxy
    integer(2) logic
    !��������� �������.
    real, parameter :: av_t = 1.0 !����� ���������� ���������� �������. ����� ��������� � �������, ��� ��� �������� ������ = 1 �.�. ����� �������� ����� �����������.
    real, parameter :: pi = 3.14159265358979
    real t, distance
    real D !����������� ��������
    !��������� ���������
    integer NumberOfParticles
    real, parameter :: MaxT = 100.0 !����� � �.�. � ������� �������� ���������� ���������.
    real, parameter :: dt = 0.1
    !��������.
    integer i, j
    !�������� ���� �������.
    type particle
        real x  !���������� .
        real y
        real z
        real(16) vx !������������ ��������.
        real(16) vy
        real(16) vz
        real(16) time2strike
    end type particle
    type(particle) part
    !���������������.
    integer numb_of_part_least_time
    integer(8) count
    type(particle), allocatable :: Particles(:)
    real l, temp, angle
    real(16) least_time_int
    real(16) sum_diff, sum_dist
    !����� ��� ��������� �������� �����.
    call RANDOM_SEED()
    !������ ������
    allocate(Particles(500))
    !��������� ������ ������.
    !��������� ��������� �������
    D = av_t * 1.0 / 3.0
    !�������� ������������ ���� � ���������� � �� ���� � ��������.
    call GraphicWindow()
    call GraphicAxes()
    !��������� ���� ��� ������ ����������� ������.
    open(1, file= "out.txt")
    !��� ��� �������� ������� �������� ����� �������, �� ����� ������� ���������� l �������� ����� ����� ����������.
    !��������� �������� Number0fParticles ������ � ������ ����������  ��������.
    
    logic = SetColor(4);
    call MoveTo_w(DBLE(0.0), DBLE(0.0), wxy)

do NumberOfParticles=1,250
    do i = 1,NumberOfParticles
        Particles(i).x = 0.0; Particles(i).y = 0.0; Particles(i).z = 0.0;
        call RANDOM_NUMBER(angle)
        angle = angle * pi * 2.0
        Particles(i).vx = cos(angle); Particles(i).vy = sin(angle); Particles(i).vz = 0.0;
        Particles(i).time2strike = Get_Random_Time2Strike()
    end do
    sum_diff = 0; count = 0
    t = 0
    do while (t < MaxT)
        sum_dist = 0
        least_time_int = MaxT
        numb_of_part_least_time = 0
        do j=1,NumberOfParticles
            if (Particles(j).time2strike < least_time_int) then
                least_time_int = Particles(j).time2strike
                numb_of_part_least_time = j
            end if
        end do
        do j=1,NumberOfParticles
            Particles(j).x = Particles(j).x + Particles(j).vx * least_time_int
            Particles(j).y = Particles(j).y + Particles(j).vy * least_time_int
            Particles(j).z = Particles(j).z + Particles(j).vz * least_time_int
            Particles(j).time2strike = Particles(j).time2strike - least_time_int
            if (j == numb_of_part_least_time) then
                Particles(j).time2strike = Get_Random_Time2Strike()
                Particles(j) = ConvertPart( Particles(j), Get_Random_CosTetha(), Get_Random_Phi() )
            end if
        end do
        do j=1,NumberOfParticles
            sum_dist = sum_dist + sqrt( (Particles(j).x)**2 + (Particles(j).y)**2 + (Particles(j).z)**2 )
        end do
        t = t + least_time_int
        if(t > MaxT / 2.0) then
            sum_diff = sum_diff + abs( sum_dist / DBLE(NumberOfParticles) - sqrt(5*D*t) )
            count = count + 1
        end if
    end do
    logic = LineTo_w( DBLE(NumberOfParticles), DBLE(sum_diff / DBLE(count) ) )
    write(1,*) NumberOfParticles, " ", sum_diff / DBLE(count)
end do

    close(1)
    
    contains
    
    real function Get_Random_Time2Strike() !��������� ���������� ������� ������� �������.
        real r
        call RANDOM_NUMBER(r)
        Get_Random_Time2Strike = -av_t * log( 1 - r )
    end function Get_Random_Time2Strike
    
    real function Get_Random_CosTetha()    !��������� ���� � ��������� ���������� ������.
        real r
        call RANDOM_NUMBER(r)
        Get_Random_CosTetha = 1 - 2 * r
    end function Get_Random_CosTetha
    
    real function Get_Random_Phi()      !���������� ���� � ��������� ����������������� ����� ���������� ������.
        real r
        call RANDOM_NUMBER(r)
        Get_Random_Phi = 2 * pi * r
    end function Get_Random_Phi
    
    type(particle) function ConvertPart( old_part, cos_tetha,  phi ) !���������� ������ ������� �������� �� ����� ����������.
        type(particle) old_part
        real cos_tetha, phi
        ConvertPart.x = old_part.x; ConvertPart.y = old_part.y; ConvertPart.z = old_part.z;
        ConvertPart.vx = cos_tetha * old_part.vx - ( old_part.vy * sin(phi) - old_part.vx * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vy = cos_tetha * old_part.vy + ( old_part.vx * sin(phi) + old_part.vy * old_part.vz * cos(phi) ) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.vz = cos_tetha * old_part.vz - ( 1 - (old_part.vz)**2 ) * cos(phi) * sqrt ( (1 - (cos_tetha)**2) / (1 - (old_part.vz)**2) )
        ConvertPart.time2strike = old_part.time2strike
    end function ConvertPart
    
    subroutine GraphicWindow() !������ ���� ��� ������ �������. ��� - �����, ������� - ������, ���� ���� - �����.
        logical(2) bool2
        integer Pxl, Pxr, Pyu, Pyd
        Pxl = 100; Pyl = 50; Pxr = 900; Pyr = 650
        bool2 = SetBkColor(15); call ClearScreen(0) !������� ����� ������ � �����.
        bool2 = SetColor(0); bool2 = Rectangle($GBORDER,Pxl-1, Pyl-1, Pxr+1, Pyr+1) !�������� ������ ����.
        call SetViewPort(Pxl, Pyl, Pxr, Pyr) !�������� ������� �������.
        bool2 = SetBkColor(15); call ClearScreen(1) !������� ������� ������� ��� ������ � �����. 1 - ������ �� ��� ������� �������.
    end subroutine GraphicWindow
    
    subroutine GraphicAxes()
        real xl, yl, xr, yr, scale_width
        real x, y
        xl = -0.1; yl = -0.1; xr = 250; yr = 5.0; scale_width = 0.1 !����������� ������ ��������� ������ ���������.
        bool2 = SetWindow(.TRUE., DBLE(xl), DBLE(yl), DBLE(xr), DBLE(yr))
        x = xl
        do while (ceiling(x) <= floor(xr)) !����������� ����� �������.
            call MoveTo_w(DBLE(ceiling(x)), DBLE(0.0 - scale_width), wxy)
            bool2 = LineTo_w(DBLE(ceiling(x)), DBLE(0.0 + scale_width))
            x = x + 1
        end do
        y = yl
        do while (ceiling(y) <= floor(yr)) !����������� ����� �������.
            call MoveTo_w(DBLE(0.0 - scale_width), DBLE(ceiling(y)), wxy)
            bool2 = LineTo_w(DBLE(0.0 + scale_width), DBLE(ceiling(y)))
            y = y + 1
        end do
        bool2 = SetColor(4) !������ ���� ���.
        call MoveTo_w(DBLE(xl), DBLE(0.0), wxy)
        bool2 = LineTo_w(DBLE(xr), DBLE(0.0))
        call MoveTo_w(DBLE(0.0), DBLE(yl), wxy)
        bool2 = LineTo_w(DBLE(0.0), DBLE(yr))
    end subroutine GraphicAxes
    
end