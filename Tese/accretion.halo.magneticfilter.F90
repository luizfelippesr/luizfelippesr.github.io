!! Copyright 2010, 2011 Luiz Felippe S. Rodrigues <felippe@astro.iag.usp.br>
!!    
!!    This file is an addition to the semi-analytical galaxy formation model
!!    Galacticus (http://sites.google.com/site/galacticusmodel/).
!!    
!!    This module is based on the "Simple" halo gas accretion module in Galacticus
!!    original code, developed by Andrew Benson.
!!
!!    This module is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!    
!!    This module is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements calculations of baryonic accretion onto halos using a Gnedin (2000) formula corrected by the influence of primordial magnetic fields (Rodrigues et al. 2010, de Souza et al. 2010).

module Accretion_Halos_MagneticFiltering
  !% Implements calculations of baryonic accretion onto halos using a Gnedin 2000 filtering mass formulation
  use FGSL
  use Radiation_Structure
  use Abundances_Structure
  use Numerical_Interpolation
  use ISO_Varying_String
  use, intrinsic :: ISO_C_Binding                             
  private 
  public :: Accretion_Halos_MagneticFiltering_Initialize, print_Mf

  ! Integration tolerance
  double precision, parameter ::  tolerance=1.0d-6

  type(varying_string) :: filteringMethod

  !% Parameters controlling reionization and magnetic fields
  double precision :: reionizationStartsRedshift,reionizationCompletetedRedshift, reionizationExponentAlpha
  double precision :: primordialMagneticFieldIntensity,reionizationTemperatureISM
  double precision :: gammaFilteringFactor

  
  ! Module global variable used in integration
  double precision :: integration_endpoint_expansion_factor
  ! Module global variable used in integration
  double precision ::  integration_small_endpoint
  
  ! Index of Solar abundance pattern.
  integer          :: abundanceIndexSolar

  ! Internal record of the number of molecules being tracked.
  integer          :: moleculesCount

  ! Zero abundance structure.
  type(abundancesStructure) :: zeroAbundances
  ! Radiation structure.
  type(radiationStructure) :: radiation
  !$omp threadprivate(radiation)


  type FilteringMassTableType
    !% Structure that holds the Filtering Mass tabulated values.
    ! number of points
    integer                                     :: nPoints 
    ! vectors, containing points
    double precision, allocatable, dimension(:) :: x, y
    ! GSL structure to help finding the index
    type(fgsl_interp_accel)                     :: interpolationAccelerator
    ! GSL structure storing relevant (whatever...) information
    type(fgsl_interp)                           :: interpolationObject
    ! Tells whether the table is being accessed for the first time or not
    logical                                     :: reset=.true.
  end type FilteringMassTableType
            
  type(FilteringMassTableType), save     :: FilteringMassTable


contains 

  !# <accretionHalosMethod>
  !#  <unitName>Accretion_Halos_MagneticFiltering_Initialize</unitName>
  !# </accretionHalosMethod>
  subroutine Accretion_Halos_MagneticFiltering_Initialize(accretionHalosMethod,Halo_Baryonic_Accretion_Rate_Get &
       &,Halo_Baryonic_Accreted_Mass_Get,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get &
       &,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accreted_Abundances_Get&
       &,Halo_Baryonic_Accretion_Rate_Molecules_Get,Halo_Baryonic_Accreted_Molecules_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    use Cosmology_Functions
    use Atomic_Data
    use Molecular_Abundances_Structure
    implicit none
    type(varying_string),          intent(in)    :: accretionHalosMethod
    procedure(),          pointer, intent(inout) :: Halo_Baryonic_Accretion_Rate_Get,Halo_Baryonic_Accreted_Mass_Get &
         &,Halo_Baryonic_Failed_Accretion_Rate_Get,Halo_Baryonic_Failed_Accreted_Mass_Get&
         &,Halo_Baryonic_Accretion_Rate_Abundances_Get,Halo_Baryonic_Accreted_Abundances_Get &
         &,Halo_Baryonic_Accretion_Rate_Molecules_Get,Halo_Baryonic_Accreted_Molecules_Get
    
    integer filteringMassTablePoints
    
    if (accretionHalosMethod == 'MagneticFiltering') then
       Halo_Baryonic_Accretion_Rate_Get            => Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get
       Halo_Baryonic_Accreted_Mass_Get             => Halo_Baryonic_Accreted_Mass_MagneticFiltering_Get
       Halo_Baryonic_Failed_Accretion_Rate_Get     => Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get
       Halo_Baryonic_Failed_Accreted_Mass_Get      => Halo_Baryonic_Failed_Accreted_Mass_MagneticFiltering_Get
       Halo_Baryonic_Accretion_Rate_Abundances_Get => Halo_Baryonic_Accretion_Rate_Abundances_MagneticFiltering_Get
       Halo_Baryonic_Accreted_Abundances_Get       => Halo_Baryonic_Accreted_Abundances_MagneticFiltering_Get
       Halo_Baryonic_Accretion_Rate_Molecules_Get  => Halo_Baryonic_Accretion_Rate_Molecules_MagneticFiltering_Get
       Halo_Baryonic_Accreted_Molecules_Get        => Halo_Baryonic_Accreted_Molecules_MagneticFiltering_Get
       ! Get the index of the solar composition abundance pattern.
       abundanceIndexSolar=Abundance_Pattern_Lookup(abundanceName="solar")
       ! Get a count of the number of molecules being tracked.
       moleculesCount=Molecules_Property_Count()
       ! Create a structure with zero abundances.
       call zeroAbundances%metallicitySet(0.0d0,adjustElements=adjustElementsReset,abundanceIndex=abundanceIndexSolar)
       ! Read parameters.
       !@ <inputParameter>
       !@   <name>gammaFilteringFactor</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Macciò's parameter to correct filtering mass. 
       !@    Use negative values to select function gamma(z)=(1+z)**1.1/11.8 instead of a constant gamma.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("gammaFilteringFactor",gammaFilteringFactor,defaultValue=1d0)
       !@ <inputParameter>
       !@   <name>filteringMethod</name>
       !@   <defaultValue>matter</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    When this parameter is set to 
       !@    - 'matter', a matter dominated Universe is assumed, allowing the use of a simplified expression to the Jeans mass.
       !@    - 'DEcut', a sharp transition from a dark matter dominated Universo to a dark energy dominated Universe is assumed (this should overestimate the impact of dark energy).
       !@    the comparison between 'DEcut' and 'matter' runs has shown only negligible differences, thus, the implementation of a complete solution for the filtering wavenumber was not completed (since it was unnecessary).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("filteringMethod",filteringMethod,defaultValue='matter')
        !@ <inputParameter>
       !@   <name>reionizationStartsRedshift</name>
       !@   <defaultValue>11.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Redshift when HII regions begin to form.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationStartsRedshift",reionizationStartsRedshift,defaultValue= 11.0d0)
       !@ <inputParameter>
       !@   <name>reionizationCompletetedRedshift</name>
       !@   <defaultValue>10.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Redshift when the universe get completely reionized.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationCompletetedRedshift",reionizationCompletetedRedshift,defaultValue=10.0d0)
       !@ <inputParameter>
       !@   <name>reionizationExponentAlpha</name>
       !@   <defaultValue>6.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Parameter that controls the rate of growth of extragalactic UV flux.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationExponentAlpha",reionizationExponentAlpha,defaultValue=6.0d0)
       !@ <inputParameter>
       !@   <name>reionizationTemperatureISM</name>
       !@   <defaultValue>1.0d4</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Temperature of the ISM during reionization epoch
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("reionizationTemperatureISM",reionizationTemperatureISM,defaultValue=1.0d4)
       !@ <inputParameter>
       !@   <name>filteringMassTablePoints</name>
       !@   <defaultValue>500</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Number of points used in the calculation of the filtering mass lookup table.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("filteringMassTablePoints",filteringMassTablePoints,defaultValue=500)
       !@ <inputParameter>
       !@   <name>primordialMagneticFieldIntensity</name>
       !@   <defaultValue>1.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@    Average primordial magnetic field intensity in nG.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter("primordialMagneticFieldIntensity",primordialMagneticFieldIntensity,defaultValue=0.0d0)
       
       ! Define the radiation structure.
       call radiation%define([radiationTypeCMB])


       call generate_Filtering_Mass_Table(filteringMassTablePoints)

    end if
    return
  end subroutine Accretion_Halos_MagneticFiltering_Initialize

  double precision function Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass, gnedinSupressionFactor
    double precision                         :: presentExp_Fact
    presentExp_Fact= Expansion_Factor(Tree_Node_Time(thisNode))     
  
    if (thisNode%isSatellite()) then
       Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get=0.0d0
    else

       gnedinSupressionFactor=(1.0d0+0.26d0*Filtering_Mass(presentExp_Fact)/Tree_Node_Mass(thisNode))**(-3.0d0) 
!       print *, 1/presentExp_Fact-1, gnedinSupressionFactor

       Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode) &
        & * gnedinSupressionFactor
       
       unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
       
       if (unaccretedMass > 0.0d0) then
             growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
             Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get=Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get+unaccretedMass &
              & *growthRate*gnedinSupressionFactor
       end if
    end if
    return
  end function Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get

  double precision function Halo_Baryonic_Accreted_Mass_MagneticFiltering_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    double precision                       :: presentExp_Fact

    presentExp_Fact= Expansion_Factor(Tree_Node_Time(thisNode))

    Halo_Baryonic_Accreted_Mass_MagneticFiltering_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)/&
        & (1.0d0+0.26d0*Filtering_Mass(presentExp_Fact)/Tree_Node_Mass(thisNode))**3.0d0
     
    return
  end function Halo_Baryonic_Accreted_Mass_MagneticFiltering_Get
  
  double precision function Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get(thisNode)
    !% Computes the baryonic accretion rate onto {\tt thisNode}.
    use Tree_Nodes
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions    
    

    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: growthRate,unaccretedMass
    double precision                         :: presentExp_Fact, gnedinSupressionFactor

    presentExp_Fact= Expansion_Factor(Tree_Node_Time(thisNode))
     
    if (thisNode%isSatellite()) then
       Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get=0.0d0
    else
       gnedinSupressionFactor=(1.0d0+0.26d0*Filtering_Mass(presentExp_Fact)/Tree_Node_Mass(thisNode))**(-3.0d0)

       Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get=(Omega_b()/Omega_0())*Tree_Node_Mass_Accretion_Rate(thisNode)&
          & *(1.0d0 -gnedinSupressionFactor)

       unaccretedMass=Tree_Node_Hot_Halo_Unaccreted_Mass(thisNode)
       if (unaccretedMass > 0.0d0) then
          growthRate=Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)
          Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get=Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get &
              & -unaccretedMass*growthRate*gnedinSupressionFactor
       end if
    end if
    return
  end function Halo_Baryonic_Failed_Accretion_Rate_MagneticFiltering_Get
 
  double precision function Halo_Baryonic_Failed_Accreted_Mass_MagneticFiltering_Get(thisNode)
    !% Computes the mass of baryons accreted into {\tt thisNode}.
    use Tree_Nodes
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions

    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    double precision                       :: presentExp_Fact

    presentExp_Fact= Expansion_Factor(Tree_Node_Time(thisNode))

    Halo_Baryonic_Failed_Accreted_Mass_MagneticFiltering_Get=(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)*(1.0d0 -&
        & 1/(1.0d0+0.26d0*Filtering_Mass(presentExp_Fact)/Tree_Node_Mass(thisNode))**3.0d0)
    return
  end function Halo_Baryonic_Failed_Accreted_Mass_MagneticFiltering_Get
  
  subroutine Halo_Baryonic_Accretion_Rate_Abundances_MagneticFiltering_Get(thisNode,accretionRateAbundances)
    !% Computes the rate of mass of abundance accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium.
    use Tree_Nodes
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(abundancesStructure), intent(inout)          :: accretionRateAbundances

    ! Assume zero metallicity.
    accretionRateAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accretion_Rate_Abundances_MagneticFiltering_Get
  
  subroutine Halo_Baryonic_Accreted_Abundances_MagneticFiltering_Get(thisNode,accretedAbundances)
    !% Computes the mass of abundances accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Tree_Nodes
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(abundancesStructure), intent(inout)          :: accretedAbundances

    ! Assume zero metallicity.
    accretedAbundances=zeroAbundances
    return
  end subroutine Halo_Baryonic_Accreted_Abundances_MagneticFiltering_Get
  
  subroutine Halo_Baryonic_Accretion_Rate_Molecules_MagneticFiltering_Get(thisNode,accretionRateMolecules)
    !% Computes the rate of mass of molecules accretion (in $M_\odot/$Gyr) onto {\tt thisNode} from the intergalactic medium. Assumes a
    !% primordial mixture of hydrogen and helium and that accreted material is in collisional ionization equilibrium at the virial
    !% temperature.
    use Tree_Nodes
    use Molecular_Abundances_Structure
    implicit none
    type(treeNode),                     intent(inout), pointer :: thisNode
    type(molecularAbundancesStructure), intent(inout)          :: accretionRateMolecules
    double precision                                           :: massAccretionRate

    ! Return immediately if no molecules are being tracked.
    if (moleculesCount == 0) return

    ! Ensure that molecules are reset to zero.
    call accretionRateMolecules%reset()

    ! Get the total mass accretion rate onto the halo.
    massAccretionRate=Halo_Baryonic_Accretion_Rate_MagneticFiltering_Get(thisNode)
    
    ! Get the mass accretion rates.
    call Get_Molecular_Masses(thisNode,massAccretionRate,accretionRateMolecules)

    return
  end subroutine Halo_Baryonic_Accretion_Rate_Molecules_MagneticFiltering_Get
  
  subroutine Halo_Baryonic_Accreted_Molecules_MagneticFiltering_Get(thisNode,accretedMolecules)
    !% Computes the mass of molecules accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Tree_Nodes
    use Molecular_Abundances_Structure
    implicit none
    type(treeNode),                     intent(inout), pointer :: thisNode
    type(molecularAbundancesStructure), intent(inout)          :: accretedMolecules
    double precision                                           :: massAccreted

    ! Return immediately if no molecules are being tracked.
    if (moleculesCount == 0) return

    ! Ensure that molecules are reset to zero.
    call accretedMolecules%reset()

    ! Total mass of material accreted.
    massAccreted=Halo_Baryonic_Accreted_Mass_MagneticFiltering_Get(thisNode)

    ! Get the masses of molecules accreted.
    call Get_Molecular_Masses(thisNode,massAccreted,accretedMolecules)

    return
  end subroutine Halo_Baryonic_Accreted_Molecules_MagneticFiltering_Get
  
  subroutine Get_Molecular_Masses(thisNode,massAccreted,molecularMasses)
    !% Compute the masses of molecules accreted (in $M_\odot$) onto {\tt thisNode} from the intergalactic medium.
    use Tree_Nodes
    use Abundances_Structure
    use Cosmological_Parameters
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Atomic
    use Numerical_Constants_Math
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Physical
    use Ionization_States
    use Molecular_Abundances_Structure
    implicit none
    type(treeNode),                     intent(inout), pointer :: thisNode
    double precision,                   intent(in)             :: massAccreted
    type(molecularAbundancesStructure), intent(out)            :: molecularMasses
    double precision                                           :: massToDensityConversion,temperature,numberDensityHydrogen&
         &,electronsDensity ,hydrogensAtomicDensity,hydrogensCationDensity
!     type(radiationStructure)                                   :: radiation
    type(molecularAbundancesStructure), save                   :: molecularDensities
    !$omp threadprivate(molecularDensities)

    ! Compute coefficient in conversion of mass to density for this node.
    massToDensityConversion=massSolar/4.0d0/Pi/(hecto*megaParsec*Dark_Matter_Halo_Virial_Radius(thisNode))**3
    
    ! Compute the temperature and density of accreting material, assuming accreted has is at the virial temperature and that the
    ! overdensity is one third of the mean overdensity of the halo.
    temperature          =Dark_Matter_Halo_Virial_Temperature(thisNode)
    numberDensityHydrogen=hydrogenByMassPrimordial*(Omega_b()/Omega_0())*Tree_Node_Mass(thisNode)*massToDensityConversion/atomicMassUnit/atomicMassHydrogen
    
    ! Set the radiation field.
    call radiation%set(thisNode)

    ! Get the molecule densities.
    call Molecular_Densities(molecularDensities,temperature,numberDensityHydrogen,zeroAbundances,radiation)

    ! Convert from densities to masses.
    call molecularDensities%numberToMass(molecularMasses)
    call molecularMasses%multiply(massAccreted*hydrogenByMassPrimordial/numberDensityHydrogen/atomicMassHydrogen)

    return
  end subroutine Get_Molecular_Masses


  subroutine generate_Filtering_Mass_Table(number_of_points)
    !% Generates a filtering mass lookup table     
    use Memory_Management
    implicit none
    integer, intent(in) :: number_of_points
    integer             :: i
    double precision    :: x
    double precision    :: step, amax, amin
    double precision, parameter :: zmax = 20
    double precision, parameter :: zmin = 0
    
    amin= 1.0d0/(1.0d0+zmax)
    amax= 1.0d0/(1.0d0+zmin) 
       
    FilteringMassTable%nPoints=number_of_points
    
    call Alloc_Array(FilteringMassTable%x,[FilteringMassTable%nPoints])
    call Alloc_Array(FilteringMassTable%y,[FilteringMassTable%nPoints])
    
    step = (amax-amin)/(FilteringMassTable%nPoints-1)
    
    do i=1,FilteringMassTable%nPoints 
      x = amin+(i-1)*step
      FilteringMassTable%x(i) = x
      FilteringMassTable%y(i) = calculate_Filtering_Mass(x)
    end do
      
    return
  end subroutine generate_Filtering_Mass_Table

  subroutine print_Mf()
    !% Prints the filtering mass (used both for debugging purposes and to produce a graph of the filtering mass)
    implicit none
    integer i,N
    double precision x, amin
    N=10000
    amin=1./(1+20.0)
    do i=1,N
      x=amin+((1.0d0-amin)/N)*i
      print *, 1/x-1, Filtering_Mass(x), Kravtsov_Filtering_Mass(x), Filtering_Mass(x)/Kravtsov_Filtering_Mass(x)
    end do
    
    call exit(1)
  end subroutine print_Mf

  double precision function Filtering_Mass(a)
    !% Finds the value of the filtering mass interpolating the lookup table
    implicit none
    double precision, intent(in)     :: a
    double precision                 :: gamma
    
!   If the user signals that he wants a variable gamma, uses Maccio's prescription
    if (gammaFilteringFactor<0) then
        gamma = a**(-1.1d0)/11.8d0
    else
        gamma = gammaFilteringFactor
    endif

    Filtering_Mass = gamma * &
                   & Interpolate(FilteringMassTable%nPoints,FilteringMassTable%x,FilteringMassTable%y,&
                   &  FilteringMassTable%interpolationObject,FilteringMassTable%interpolationAccelerator,&
                   &  a, reset=FilteringMassTable%reset)
    return
  end function Filtering_Mass
  
  double precision function calculate_Filtering_Mass(a)
    !% Computes the filtering mass
    use Numerical_Integration
    use Galacticus_Error
    use Cosmological_Parameters
    implicit none
    double precision, intent(in)     :: a
    double precision                 :: a_actual
    double precision                 :: integrTwoThirds
    type(c_ptr)                      :: parameterPointer 
    type(fgsl_function)              :: integrandFunction
    type(fgsl_integration_workspace) :: integrationWorkspace
    double precision                 :: a_DE_Dominated_Phase

!   A dark matter dominated universe is assumed 
    if (filteringMethod == 'matter') then
        a_actual=a
        
!   A universe with a sharp transition between DM dominated and DE dominated is assumed        
    else if (filteringMethod == 'DEcut') then

        a_DE_Dominated_Phase = (Omega_0()/Omega_DE())**(1d0/3d0)

!       The integral is truncated during Dark Energy domination
        a_actual = min(a,a_DE_Dominated_Phase)

!       This part was removed after the discovery that DE played a small role (comparing the results from 'DEcut' runs and 'matter' runs).  
!       else if (filteringMethod == 'full') then
!         calculate_Filtering_Mass=Filtering_Mass_alt(a)
!         return  

!   Kravtsov's (2004) fit (no magnetic fields)
    else if (filteringMethod == 'Kravtsov') 
        calculate_Filtering_Mass = Kravtsov_Filtering_Mass(a)
        return
        
    end if

    ! Sets the endpoint of the integral
    integration_endpoint_expansion_factor = a_actual

    integrTwoThirds = Integrate(0.0d0,a_actual,Filtering_Mass_Integrand_two_thirds, &
                      & parameterPointer,integrandFunction,integrationWorkspace,hasSingularities=.true., toleranceAbsolute=tolerance, toleranceRelative=tolerance)
        
    call Integrate_Done(integrandFunction,integrationWorkspace)

    calculate_Filtering_Mass=(3.0d0*integrTwoThirds/a)**(3.0d0/2.0d0)
    ! In this case, the scale factor is the normal one (not a_actual), since it acts as part of the conversion of k into \lambda
    
    return 
  end function calculate_Filtering_Mass
  

  function Filtering_Mass_Integrand_two_thirds(ef, parameterPointer) bind(c)
    implicit none
    real(c_double)         :: Filtering_Mass_Integrand_two_thirds
    real(c_double), value  :: ef
    type(c_ptr),    value  :: parameterPointer
    
    Filtering_Mass_Integrand_two_thirds=((Jeans_Mass(ef))**(2.0d0/3.0d0)) * (1.0d0 - dsqrt(ef/integration_endpoint_expansion_factor) )
    
    return
  end function Filtering_Mass_Integrand_two_thirds

  double precision function Kravtsov_Filtering_Mass(exp_factor)
    !% Computes the filtering mass using the analytical fit from Kravtsov (2004)
    use Cosmological_Parameters
    implicit none
    double precision,   intent(in)    :: exp_factor 
    double precision                  :: MJ0, mu, f,al, a0, ar, hzinho, a

    hzinho=Little_H_0()
    mu = 0.59d0
    al = reionizationExponentAlpha
    a = exp_factor
    a0 = 1.0d0/(1.0d0+reionizationStartsRedshift)
    ar = 1.0d0/(1.0d0+reionizationCompletetedRedshift)

    MJ0 = 2.5d11/(hzinho  *omega_0()**(1.0d0/2.0d0)*mu**(3.0d0/2.0d0))

    if ( a <= a0 ) then
       f=3.0d0 * a**(al+1.0d0) / ( (2.0d0+al) * (5.0d0+2.0d0*al) * a0**al)
    elseif ( a <= ar ) then
       f=(3.0d0/a)*( a0**2.0d0 * (1.0d0/(2.0d0+al)-2.0d0*(a0/a)**(1.0d0/2.0d0)/(5.0d0+2.0d0*al))+(a**2.0d0)/10.0d0&
            & - (a0**2.0d0)/10.0d0 * (5.0d0 - 4.0d0*(a0/a)**(1.0d0/2.0d0)))
    else
       f=(3.0d0/a) * ( a0**2.0d0 * (1.0d0/(2.0d0+al)-2.0d0*(a0/a)**(1.0d0/2.0d0)/(5.0d0+2.0d0*al)) &
            & - (a0**2.0d0)/10.0d0 * (5.0d0 - 4.0d0*(a0/a)**(1.0d0/2.0d0)) &
            & + (ar**2.0d0)/10.0d0 * (5.0d0 - 4.0d0*(ar/a)**(1.0d0/2.0d0)) &
            & + a*ar/3.0d0 - (ar**2)/3 * (3.0d0-2.0d0*(ar/a)**(1.0d0/2.0d0)))
    endif
    Kravtsov_Filtering_Mass=MJ0*f**(3.0d0/2.0d0)
    return
  end function Kravtsov_Filtering_Mass

  double precision function Jeans_Mass(exp_factor)
  !% Computes the Jeans mass that will be later used for the calculation of the filtering mass
    use Cosmological_Parameters
    use Numerical_Constants_Physical
    use Numerical_Constants_Atomic
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in) :: exp_factor
    double precision             :: Coefficient, thermalTerm, magneticTerm
    double precision             :: mu, B2
!     double precision, parameter  :: convert_nG2_to_kg_per_m_s2=nano**2.0*kilo**(-1.0)/centi
    double precision, parameter  :: convert_nG2_to_kg_per_m_s2=1d-19
!     double precision, parameter  :: convert_m2_s2_to_km2_s2=kilo**(-2.0)
    double precision, parameter  :: convert_m2_s2_to_km2_s2=1d-6

    mu = 0.59d0  ! Mean molecular weight of fully ionized gas 
    B2 = primordialMagneticFieldIntensity*primordialMagneticFieldIntensity * convert_nG2_to_kg_per_m_s2
    
    magneticTerm = convert_m2_s2_to_km2_s2* B2 * 2d0/3d0 * gravitationalConstant /(Omega_0()*(H_0_invGyr()/gigaYear)**2d0)

    Coefficient = 8d0/3d0*(Pi**3)*dsqrt(2.0d0/3d0)*dsqrt(1/Omega_0()) / (H_0()*gravitationalConstantGalacticus)

    thermalTerm = convert_m2_s2_to_km2_s2* (5.0d0/3.0d0) * boltzmannsConstant*exp_factor*ISM_Temperature(exp_factor)/(mu*massHydrogenAtom)
    
    Jeans_Mass = Coefficient * (thermalTerm + magneticTerm)**(3.0d0/2.0d0)    
    return

  end function Jeans_Mass

  double precision function ISM_Temperature(a)
  !% Parametrization of the temperature used in the Jeans Mass calculation (using the same parametrization as Kravtsov
    implicit none
    double precision, intent(in)   :: a
    double precision               :: T4, as, al,ar
    al = reionizationExponentAlpha
    as = 1.0d0/(1.0d0+reionizationStartsRedshift)
    ar = 1.0d0/(1.0d0+reionizationCompletetedRedshift)
    if ( a<= as ) then
      T4 = (a/as)**reionizationExponentAlpha
    elseif ( a<=ar .and. a>as ) then
      T4 = 1.0d0
    elseif (a>ar ) then
      T4 = (ar/a)
    else
      call exit(1)
    endif
    ISM_Temperature = reionizationTemperatureISM*T4
    return
  end function ISM_Temperature



! ------------------------------------------------------------------------------------
!  Below this point there is a tentative implementation of a full resolution of the 
!  differential equations for the filtering wavelength. Since no big differences in the
!  luminosity and mass functions could be seen when comparing the 'DEcut'-case (sharp
!  transition from a fully DM dominated universe to a DE dominated one) and 
!  'matter'-case (DM domination), the implementation of the detailed case was 
!  interrupted (and the its many bugs and numerical mistakes were NOT corrected).
!  -------------------------------------------------------------------------------------

! double precision function integral_small(x)
!     !% Computes the filtering mass
!     use Numerical_Integration
!     use Galacticus_Error
! 
!     implicit none
!     double precision, intent(in)     :: x
!     type(c_ptr)                      :: parameterPointer
!     type(fgsl_function)              :: integrandFunction
!     type(fgsl_integration_workspace) :: integrationWorkspace
! 
! 
!     integral_small = Integrate(x,integration_small_endpoint,integrand_small, &
!                       & parameterPointer,integrandFunction,integrationWorkspace,hasSingularities=.true., toleranceAbsolute=tolerance, toleranceRelative=tolerance)
! 
!     call Integrate_Done(integrandFunction,integrationWorkspace)
! 
!     return
!   end function integral_small
! 
! 
!   function integrand_small(y, parameterPointer) bind(c)
!     use Cosmology_Functions
!     implicit none
!     real(c_double)         :: integrand_small
!     real(c_double), value  :: y
!     type(c_ptr),    value  :: parameterPointer
! 
! 
!     integrand_small=(Expansion_Factor(y))**(-2d0)
! 
!     return
!   end function integrand_small
! 
!   double precision function integral_big(t)
!     !% Computes the filtering mass
!     use Numerical_Integration
!     use Galacticus_Error
! 
!     implicit none
!     double precision, intent(in)     :: t
!     type(c_ptr)                      :: parameterPointer
!     type(fgsl_function)              :: integrandFunction
!     type(fgsl_integration_workspace) :: integrationWorkspace
!     double precision, parameter      :: smallestTime=5d-4
!     integration_small_endpoint=t
! 
!     
!     integral_big = Integrate(smallestTime,t,integrand_big, &
!                       & parameterPointer,integrandFunction,integrationWorkspace,hasSingularities=.true., toleranceAbsolute=tolerance, toleranceRelative=tolerance)
! 
!     call Integrate_Done(integrandFunction,integrationWorkspace)
! 
!     return
!   end function integral_big
! 
!   double precision function Ddot(t)
!     use Linear_Growth
!     use Cosmology_Functions
!     implicit none
!     double precision, intent(in)   :: t
! 
!     Ddot = Expansion_Rate(Expansion_Factor(t)) &
!        &  * Linear_Growth_Factor_Logarithmic_Derivative(time=t) &
!        &  * Linear_Growth_Factor(time=t)
!     return
!   end function Ddot
!   
!   function integrand_big(x, parameterPointer) bind(c)
!     use Linear_Growth
!     use Cosmology_Functions
!     implicit none
!     real(c_double)         :: integrand_big
!     real(c_double), value  :: x
!     type(c_ptr),    value  :: parameterPointer
!     double precision       :: D2dot
!     double precision       :: h
! 
! !     h=max(1d-5,1d-8*x)
!     h=1d-12*x
! !   Second derivative with respect to time
!     D2dot= ( Ddot(x+h) - Ddot(x-h) )/(2d0*h)
! 
!     print *, 1/Expansion_Factor(x) -1, x,Ddot(x+h)/Ddot(x-h)
! 
!     if (abs(D2dot)<1d-10) D2dot=0
!     if (D2dot<0 .or. Expansion_Rate(Expansion_Factor(x))<0 .or. Ddot(x)<0 ) then
! !         h=1d-7*x
!         print *, 'z=',1/Expansion_Factor(x) -1, 'ddot',Ddot(x-h),Ddot(x+h), 'd2dot', D2dot, 'H',Expansion_Rate(Expansion_Factor(x))
!         stop
!     endif
!     
!     integrand_big=(Expansion_Factor(x))**(2d0) * ( D2dot + &
!                    & 2d0*Expansion_Rate(Expansion_Factor(x)) * Ddot(x)) &
!                    & /k2_Jeans(x) !* integral_small(x)
! 
!     return
!   end function integrand_big
! 
!   double precision function rho(exp_factor)
!     use Cosmological_Parameters
!     use Cosmology_Functions
!     use Numerical_Constants_Math
!     use Numerical_Constants_Physical
!     implicit none
!     double precision, intent(in)   :: exp_factor
! 
!     rho=omega_0()*H_0_invGyr()**2d0 / (8d0*Pi*gravitationalConstantGalacticus) * exp_factor**(-3d0)
! 
!   end function rho
! 
! 
!   double precision function k2_Jeans(x)
!     use Ideal_Gases_Thermodynamics
!     use Numerical_Constants_Math
!     use Cosmology_Functions
!     implicit none
!     double precision, intent(in)   :: x
!     double precision               :: sound_speed_squared, alfven_speed_squared
!     double precision               :: exp_factor
! 
!     exp_factor = Expansion_Factor(x)
! 
!     sound_speed_squared = (Ideal_Gas_Sound_Speed(ISM_Temperature(exp_factor)))**2d0
! 
!     alfven_speed_squared = 0d0
! 
!     k2_Jeans= 4d0*Pi*rho(exp_factor)/(alfven_speed_squared + sound_speed_squared)
! print *, 'Jeans',1d0/dsqrt(k2_Jeans), 'z', 1/exp_factor-1
!   end function k2_Jeans
! 
! 
!   double precision function Filtering_Mass_alt(exp_factor)
!     use Linear_Growth
!     use Cosmological_Parameters
!     use Cosmology_Functions
!     use Numerical_Constants_Math
!     implicit none
!     double precision, intent(in)     :: exp_factor
!     double precision                 :: Filtering_Length, t
!     double precision IB
! !   tempo em Ganos
!     t=Cosmology_Age(exp_factor)
! 
! !     Galacticus' M_Solar, Mpc, km/s unit system
! 
! !   M_solar/Mpc^3
! 
!      IB = integral_big(t)
!      print *, IB
!      
!     Filtering_Length = 2d0*Pi*exp_factor* dsqrt(Linear_Growth_Factor(time=t)/IB)
! 
!     Filtering_Mass_alt = 4d0/3d0 * Pi *rho(exp_factor) * (Filtering_Length)**(3d0)
! 
!   end function Filtering_Mass_alt
  
end module Accretion_Halos_MagneticFiltering
