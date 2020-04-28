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
    real, parameter :: MaxT = 1000.0 !����� � �.�. � ������� �������� ���������� ���������.
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
    end type particle
    type(particle) part
    !���������������.
    real l, temp, angle
    real(16) sum_diff
    real(16) SUM_DIFF_K
    real dt
    !����� ��� ��������� �������� �����.
    call RANDOM_SEED()
    !��������� ��������� �������
    D = av_t * 1.0 / 3.0    !!!��� ��� �������� ������� �������� ����� �������, �� ����� ������� ���������� l �������� ����� ����� ����������!!!
    !�������� ������������ ���� � ���������� � �� ���� � ��������.
    call GraphicWindow()
    call GraphicAxes()
    !��������� ���� ��� ������ ����������� ������.
    open(1, file= "out.txt")
    
logic = SetColor(4);
call MoveTo_w(DBLE(0.0), DBLE(0.0), wxy)
do NumberOfParticles=1,300
    call RANDOM_SEED()
    !������� ���������� ������ �� ������ ��� ������ ������� �� K ������.
    SUM_DIFF_K = 0
    do i=1,NumberOfParticles
        sum_diff = 0
        !������� ���������, ������� i-�� ������� ������ �� ����� Tmax.
        t = 0
        !��������� ��������� ���� �������.
        part.x = 0.0; part.y = 0.0; part.z = 0.0;
        part.vx = 1.0; part.vy = 0.0; part.vz = 0.0;
        !������� ���������, �� ������� �������� �� ������ ��������� i-�� ������� �� ����� MaxT.
        do while (t < MaxT)
            dt = Get_Random_Time2Strike()
            if (t + dt < MaxT) then
                part.x = part.x + part.vx * dt
                part.y = part.y + part.vy * dt
                part.z = part.z + part.vz * dt
                part = ConvertPart( part, Get_Random_CosTetha(), Get_Random_Phi() )
            else
                dt = (MaxT - t)
                part.x = part.x + part.vx * dt
                part.y = part.y + part.vy * dt
                part.z = part.z + part.vz * dt
            end if
            t = t + dt
        end do
        sum_diff = sum_diff + abs( sqrt( (part.x)**2 + (part.y)**2 + (part.z)**2 ))
    end do
    !������� �������.
    sum_diff = abs(sum_diff / DBLE(NumberOfParticles) - DBLE(sqrt(6*D*MaxT)))
    !������ ��������������� ����� �� �������.
    logic = LineTo_w( DBLE(NumberOfParticles), DBLE(sum_diff) )
    write(1,*) (DBLE(NumberOfParticles)), " ", sum_diff
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
        xl = -0.1; yl = -0.1; xr = 300; yr = 100.0; scale_width = 0.1 !����������� ������ ��������� ������ ���������.
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