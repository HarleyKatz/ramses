!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! EMISSION LINE LUMINOSITY FUNCTIONS !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE movie_lines_module
    use amr_commons, only: dp
    implicit none

    private
    public :: get_halpha_lum, get_hbeta_lum, get_OIII_5007_lum, &
            & get_OIII_4363_lum, get_OIII_4959_lum, get_OII_3726_lum, &
            & get_OII_3728_lum, get_NII_6583_lum

    real(dp):: clight  = 2.99792458d10          ![cm/s] light speed
    real(dp):: planck  = 6.626070040d-27        ![erg s] Planck's constant

CONTAINS

FUNCTION get_halpha_lum(T,ne,nHII,nHI,vol) result(lum)
    !
    ! Returns the Halpha luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nHII --> ionized hydrogen density in cm^-3
    !nHI  --> neutral hydrogen density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nHII,nHI,vol
    real(dp)::lum
    real(dp)::ener
    real(dp)::rec_emis, col_emis
    real(dp)::Z,a,b,c,d,et

    ! Recombination H alpha emissivity erg cm^3 s^-1 from Pequignot+ 1991
    ! Recombination formula from 
    ! http://articles.adsabs.harvard.edu//full/1991A%26A...251..680P/0000684.000.html
    ener = (planck * clight) / (6562.8d0 * 1d-8) 
    Z = 1.d0
    a = 2.708d0
    b = -0.648d0
    c = 1.315d0
    d = 0.523d0
    et = (1d-4) * T / (Z**2)
    rec_emis = ener * 1d-13 * Z * (a * et**b) / (1.d0 + c * et**d)
    rec_emis = rec_emis * ne * nHII * vol

    ! Collisional emission calculated from harley's fitting formula
    a = 5.01d-19
    b = 8.13d4
    c = 0.230d0
    d = 0.938d0
    col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
    col_emis = col_emis * ne * nHI * vol

    lum = rec_emis + col_emis ! erg/s

END FUNCTION get_halpha_lum
        
FUNCTION get_hbeta_lum(T,ne,nHII,nHI,vol) result(lum)
    ! Returns the Hbeta luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nHII --> ionized hydrogen density in cm^-3
    !nHI  --> neutral hydrogen density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nHII,nHI,vol
    real(dp)::lum
    real(dp)::ener
    real(dp)::rec_emis, col_emis
    real(dp)::Z,a,b,c,d,et

    ! Recombination H alpha emissivity erg cm^3 s^-1 from Pequignot+ 1991
    ! Recombination formula from 
    ! http://articles.adsabs.harvard.edu//full/1991A%26A...251..680P/0000684.000.html
    ener = (planck * clight) / (4861.4d0 * 1d-8) 
    Z = 1.d0
    a = 0.668d0
    b = -0.507d0
    c = 1.221d0
    d = 0.653d0
    et = (1d-4) * T / (Z**2)
    rec_emis = ener * 1d-13 * Z * (a * et**b) / (1.d0 + c * et**d)
    rec_emis = rec_emis * ne * nHII * vol

    ! Collisional emission calculated from harley's fitting formula
    a = 1.81d-19
    b = 9.87d4
    c = 0.237d0
    d = 0.954d0
    col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
    col_emis = col_emis * ne * nHI * vol

    lum = rec_emis + col_emis ! erg/s

END FUNCTION get_hbeta_lum
        
FUNCTION get_OIII_5007_lum(T,ne,nO3,vol) result(lum)
    !
    ! Returns the OIII_5007 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nO3  --> OIII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nO3,vol
    real(dp)::lum
    real(dp)::a,b,c,d,Tmin,Tmax,col_emis

    ! Collisional emissivity
    ! Based on Harley's fits to pyneb data
    ! This is empirically good up to n_e = 1e5 
    Tmin = 1325.7113655901096d0
    Tmax = 10000000.0d0
    lum = 0.d0
    if (T.ge.Tmin.and.T.le.Tmax) then 
        a = 8.54087422d-18
        b = 2.36679211d+04
        c = 5.03985699d-01
        d = 9.69115359d-01
        col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
        lum = col_emis * ne * nO3 * vol
    end if
END FUNCTION get_OIII_5007_lum
        
FUNCTION get_OIII_4363_lum(T,ne,nO3,vol) result(lum)
    ! Returns the OIII_4363 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nO3  --> OIII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nO3,vol
    real(dp)::lum
    real(dp)::a,b,c,d,Tmin,Tmax,col_emis

    ! Collisional emissivity
    ! Based on Harley's fits to pyneb data
    ! We fit for the ratio of 4363/5007
    lum = get_OIII_5007_lum(T,ne,nO3,vol) 

    Tmin = 2682.6957952797275d0
    Tmax = 10000000.0d0
    lum = 0.d0
    if (T.ge.Tmin.and.T.le.Tmax) then 
        a = 1.58694188d-01
        b = 4.30073857d+04
        c = -4.70853221d-04
        d = 1.03287240d+00
        col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
        lum = col_emis * lum
    end if
END FUNCTION get_OIII_4363_lum
        
FUNCTION get_OIII_4959_lum(T,ne,nO3,vol) result(lum)
    !
    ! Returns the OIII_4959 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nO3  --> OIII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nO3,vol
    real(dp)::lum

    ! Collisional emissivity
    ! We pin the result to the 5007 fitting function
    lum = get_OIII_5007_lum(T,ne,nO3,vol) / 2.98d0

END FUNCTION get_OIII_4959_lum

FUNCTION softplus(x,a,b,c)
    implicit none
    real(dp)::x,a,b,c
    real(dp)::softplus
    softplus = LOG(a + EXP(b*x - c))
END FUNCTION softplus
        
FUNCTION softplus_neg(x,a,b,c) result(fun)
    implicit none
    real(dp), intent(in)::x,a,b,c
    real(dp)::fun
    fun = a - LOG(1.0 + EXP(b*x - c))
END FUNCTION softplus_neg
        
FUNCTION get_OII_3726_lum(T,ne,nO2,vol) result(lum)
    !
    ! Returns the OII_3726 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nO3  --> OII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nO2,vol
    real(dp)::lum
    real(dp)::a,b,c,d,Tmin,Tmax,logne,col_emis

    ! Collisional emissivity
    ! Based on Harley's fits to pyneb data
    ! This is empirically good up to n_e = 1e5 
    Tmin = 1325.7113655901096d0
    Tmax = 10000000.0d0
    logne = LOG10(ne)
    lum = 0.d0
    if (T.ge.Tmin.and.T.le.Tmax) then 
        a = 10.0d0**( softplus_neg(logne,-17.0923918d0,2.63282097d0,  10.44740535d0) )   
        b = 30000.0d0 * softplus(logne,2.93457396d0,1.54750724d0, 5.90262693d0) 
        c = 10.0d0**( softplus_neg(logne,-0.30407987d0,1.821353d0,  9.44135231d0) )
        d = softplus(logne,2.6409534d0,0.92271611d0,5.83020548d0)
        col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
        lum = col_emis * ne * nO2 * vol
    end if
END FUNCTION get_OII_3726_lum
        
FUNCTION get_OII_3728_lum(T,ne,nO2,vol) result(lum)
    !
    ! Returns the OII_3728 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nO3  --> OII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nO2,vol
    real(dp)::lum
    real(dp)::a,b,c,d,Tmin,Tmax,logne,col_emis

    ! Collisional emissivity
    ! Based on Harley's fits to pyneb data
    ! This is empirically good up to n_e = 1e5 
    Tmin = 1325.7113655901096d0
    Tmax = 10000000.0d0
    logne = LOG10(ne)
    lum = 0.d0
    if (T.ge.Tmin.and.T.le.Tmax) then 
        a = 10.0d0**( softplus_neg(logne,-16.91037079d0,2.01549352d0,6.24261369d0) )   
        b = 30000.0d0 * softplus(logne,2.92460964d0,2.115382d0,7.5600845d0) 
        c = 10.0d0**( softplus_neg(logne,-0.28833215d0,1.18957554d0,5.68403389d0) )
        d = softplus(logne,2.63855695d0,1.06518539d0,6.06570019d0)
        col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
        lum = col_emis * ne * nO2 * vol
    end if
END FUNCTION get_OII_3728_lum
        
FUNCTION get_NII_6583_lum(T,ne,nN2,vol) result(lum)
    !
    ! Returns the NII_6583 luminosity of a cell in erg/s
    !
    !T    --> Temperature in K
    !ne   --> electron density in cm^-3
    !nN2  --> NII density in cm^-3
    !vol  --> cell volume in cm^3
    implicit none
    real(dp), intent(in)::T,ne,nN2,vol
    real(dp)::lum
    real(dp)::a,b,c,d,Tmin,Tmax,col_emis

    ! Collisional emissivity
    ! Based on Harley's fits to pyneb data
    ! This is empirically good up to n_e = 1e5 
    Tmin = 1048.1131341546852d0
    Tmax = 10000000.0d0
    lum = 0.d0
    if (T.ge.Tmin.and.T.le.Tmax) then 
        a = 5.93148982e-18    
        b = 1.98682195e+04 
        c = 4.92946876e-01 
        d = 9.83930826e-01
        col_emis = (a / T**c) * EXP(-1.d0 * b / T**d)
        lum = col_emis * ne * nN2 * vol
    end if
END FUNCTION get_NII_6583_lum

end module movie_lines_module