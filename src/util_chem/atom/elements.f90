!--------------------------------------------------------------------------------------------------
! MODULE: element_mod
!> @author Marin Sapunar, Ruđer Bošković Institute
!> @date October, 2017
!
!> @brief Contains lists of element data and subroutines for accessing it.
!--------------------------------------------------------------------------------------------------
module element_mod
    ! Import variables
    use global_defs
    ! Import subroutines
    use string_mod, only : tolower
    implicit none
   
    private
    public :: element_symbol
    public :: element_arstand
    public :: element_z2s
    public :: element_z2m
    public :: element_s2z
    public :: element_s2m

    character(len=2), parameter, dimension(1:103) :: element_symbol = [ &
      "H ", "He", & ! Row 1
      "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", & ! Row 2
      "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", & ! Row 3
      "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", & ! Row 4
      "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", & ! Row 4 (cont.)
      "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", & ! Row 5
      "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", & ! Row 5 (cont.)
      "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", & ! Row 6
      "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", & ! Row 6 (cont.)
      "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", & ! Row 6 (cont.)
      "Pb", "Bi", "Po", "At", "Rn",                         & ! Row 6 (cont.)
      "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", & ! Row 7
      "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"        & ! Row 7 (cont.)
      ]

    real(dp), parameter, dimension(1:103) :: element_arstand = [ &
      1.00794_dp, 4.002602_dp, 6.941_dp, 9.01218_dp, 10.811_dp, 12.011_dp, 14.00674_dp, 15.9994_dp, 18.998403_dp, 20.1797_dp, 22.989768_dp, 24.305_dp, 26.981539_dp, 28.0855_dp, 30.973762_dp, 32.066_dp, 35.4527_dp, 39.948_dp, 39.0983_dp, 40.078_dp, 44.95591_dp, 47.88_dp, 50.9415_dp, 51.9961_dp, 54.93805_dp, 55.847_dp, 58.9332_dp, 58.6934_dp, 63.546_dp, 65.39_dp, 69.723_dp, 72.61_dp, 74.92159_dp, 78.96_dp, 79.904_dp, 83.8_dp, 85.4678_dp, 87.62_dp, 88.90585_dp, 91.224_dp, 92.90638_dp, 95.94_dp, 97.9072_dp, 101.07_dp, 102.9055_dp, 106.42_dp, 107.8682_dp, 112.411_dp, 114.818_dp, 118.71_dp, 121.76_dp, 127.6_dp, 126.90447_dp, 131.29_dp, 132.90543_dp, 137.327_dp, 138.9055_dp, 140.115_dp, 140.90765_dp, 144.24_dp, 144.9127_dp, 150.36_dp, 151.965_dp, 157.25_dp, 158.92534_dp, 162.5_dp, 164.93032_dp, 167.26_dp, 168.93421_dp, 173.04_dp, 174.967_dp, 178.49_dp, 180.9479_dp, 183.84_dp, 186.207_dp, 190.23_dp, 192.22_dp, 195.08_dp, 196.96654_dp, 200.59_dp, 204.3833_dp, 207.2_dp, 208.98037_dp, 208.9824_dp, 209.9871_dp, 222.0176_dp, 223.0197_dp, 226.0254_dp, 227.0278_dp, 232.0381_dp, 231.03588_dp, 238.0289_dp, 237.048_dp, 244.0642_dp, 243.0614_dp, 247.0703_dp, 247.0703_dp, 251.0796_dp, 252.083_dp, 257.0951_dp, 258.1_dp, 259.1009_dp, 262.11_dp &
      ]


contains

    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE Element_z2s
    !> @brief Convert atomic number to element symbol.
    !----------------------------------------------------------------------------------------------
    elemental function element_z2s(z) result (sym)
        integer, intent(in) :: z
        character(len=2) :: sym

        if (z > size(element_symbol)) then
            sym = "??"
        else
            sym = element_symbol(z)
        end if
    end function element_z2s


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE Element_z2m
    !> @brief Convert atomic number to standard atomic weight.
    !----------------------------------------------------------------------------------------------
    elemental function element_z2m(z) result (m)
        integer, intent(in) :: z
        real(dp) :: m

        if (z > size(element_arstand)) then
            m = -1.0_dp
        else
            m = element_arstand(z)
        end if
    end function element_z2m


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE Element_s2z
    !> @brief Convert element symbol to atomic number.
    !----------------------------------------------------------------------------------------------
    elemental function element_s2z(sym) result (z)
        character(len=2), intent(in) :: sym
        integer :: z
        integer :: i

        do i = 1, size(element_symbol)
            if (trim(tolower(element_symbol(i))) /= trim(tolower(adjustl(sym)))) cycle
            z = i
            return
        end do
        z = -1
    end function element_s2z


    !----------------------------------------------------------------------------------------------
    ! SUBROUTINE Element_s2m
    !> @brief Convert element symbol to standard atomic weight.
    !----------------------------------------------------------------------------------------------
    elemental function element_s2m(sym) result (m)
        character(len=2), intent(in) :: sym
        real(dp) :: m

        m = element_z2m(element_s2z(sym))
    end function element_s2m


end module element_mod
