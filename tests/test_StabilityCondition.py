import pytest
import os
import sys

from sympy import symbols, expand

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.DerivedCategory.ChernCharacter import ChernCharacter
from src.DerivedCategory.StabilityCondition import HNFactor, HarderNarasimhanError, HarderNarasimhanFiltration, StabilityCondition
from src.DerivedCategory.GeometryContext import GeometryContext, DivisorData
from src.DerivedCategory.CoherentSheaf import LineBundle
from src.DerivedCategory.DerivedCategoryObject import GradedCoproductObject
from src.DerivedCategory.SphericalTwist import SphericalTwistComposition



def test_StabilityCondition_init():
    H, C, D = symbols("H C D")
    basis3 = [H, C, D]
    tensor_data = {
        (H, H): 4,
        (C, C): -1,
        (D, D): -1,
        (H, C): 1,
        (C, D) : 0,
        (D, H) : 2
    }

    divisor_data = DivisorData(basis=basis3, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data)

    sigma1 = StabilityCondition(geometry_context, 1,2,3,4)
    
    tensor_data_2 = {
        (H,): 1
    }
    divisor_data_2 = DivisorData(basis=[H], top_intersection_form=tensor_data_2)
    geometry_context2 = GeometryContext(catagory='P1', divisor_data=divisor_data_2)

    sigma2 = StabilityCondition(geometry_context2, complex(5,10))

    tensor_data_3 = {
        (H,H) : 1
    }
    divisor_data_3 = DivisorData(basis=[H], top_intersection_form=tensor_data_3)
    geometry_context3 = GeometryContext(catagory='P2', divisor_data=divisor_data_3)
    sigma3 = StabilityCondition(geometry_context3, 15.5, 25)

    
    with pytest.raises(ValueError):
        StabilityCondition(geometry_context, 0.5, 1.5, 2.4)
    with pytest.raises(TypeError):
        StabilityCondition(geometry_context, "0.5", "2.4", "3.4", "4")
    with pytest.raises(TypeError):
        StabilityCondition(geometry_context2, 0.5)


def test_StabilityCondition_central_charge():


    #####################
    # Test K3 surface picard rank 1
    #####################

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H): 4,
    }
    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)
    sigma1 = StabilityCondition(geometry_context, 2, 3)
     
    lb1 = LineBundle(2*H, geometry_context=geometry_context)
    lb2 = LineBundle(3*H, geometry_context=geometry_context)

    ## Mukai vector of lb2 should be <1, 2H, 2H**2 + 1> 
    ## Central charge computation should be 
    ##
    ## Z_2,2(2H) = <1, 2H + 3i H, (4 - 9)*H**2/2 + 6i H**2> * <1, 2H, 2H**2 + 1>
    ##           = 4H**2 + 6i H**2 + 5/2 H**2 - 6i H**2 - 2H**2 -1
    ##           = 9/2*H**2 - 1
    ##           = 17

    assert sigma1.centralCharge(lb1) == 17

    ## Mukai vector of lb3 should be <1, 3H, 9H**2/2 + 1>
    ## Central charge computation should be
    ##
    ## Z_2,2(3H) = <1, 2H + 3i H, (4 - 9)*H**2/2 + 6i H**2> * <1, 3H, 9H**2/2 + 1>
    ##           = (6*H**2 + 9i*H**2 + 5/2*H**2 - 6i*H**2 - 9/2*H**2) - 1
    ##           = (4*H**2 + 3i*H**2) - 1
    ##           = 15 + 12i


    assert sigma1.centralCharge(lb2) == 15 + 12j


    #####################
    # Test K3 surface of picard rank 2
    #####################
    H, C = symbols("H C")
    basis = [H, C]
    tensor_data = {
        (H, H): 4,
        (C, C): -1,
        (H, C): 2,
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)
    sigma = StabilityCondition(geometry_context, 1, 1, 3)

    lb1 = LineBundle(H + 2*C, geometry_context=geometry_context)
    lb2 = LineBundle(H + 3*C, geometry_context=geometry_context)

    ## Mukai vector of lb1 should be <1, H + 2C, 0.5*H**2 + 2*C**2 + 2*H*C + 1>
    ## Central charge computation should be
    ##
    ## Z_(H + C, 3H) = < 1, H + C + 3iH, 0.5*(-8H^2 +2HC + C^2) + 3i(H^2 + HC)> * <1, H + 2C, 0.5*H**2 + 2*C**2 + 2*H*C + 1>
    ##               = (H^2 + 3HC + 2C^2 + 3i (H^2 + 2HC)) 
    ##                 - (0.5H^2 - 2C^2 - 2HC - 1)
    ##                 + 4H^2 - HC - 0.5 C^2 - 3i H^2 - 3i HC
    ## 
    ##               = 4.5 H^2 - 0.5 C^2 + 3i HC - 1

    assert sigma.centralCharge(lb1) == 17.5 + 6j


    ## Mukai vector of lb2 should be <1, H + 3C, 0.5*H**2 + 4.5*C**2 + 3*H*C + 1>
    ## Central charge computation should be

    ## Z_(H + C, 3H) = < 1, H + C + 3iH, 0.5*(-8H^2 +2HC + C^2) + 3i(H^2 + HC)> * <1, H + 3C, 0.5*H**2 + 4.5*C**2 + 3*H*C + 1>
    ##               =  H^2 + 4CH + 3C^2 + 3iH^2 + 9iHC
    ##                 - 0.5H^2 - 4.5 C**2 - 3HC - 1
    ##                  +4 H^2 - HC - 0.5C^2 - 3i H^2 - 3iHC
    ##
    ##               =  4.5 H^2 - 2C^2 + i(6HC) - 1

    assert sigma.centralCharge(lb2) == 19+12j


def test_StabilityCondition_phase():

    #####################
    # Test K3 surface picard rank 1
    #####################

    H = symbols("H")
    basis = [H]
    tensor_data = {
        (H, H): 4,
    }
    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)
    sigma1 = StabilityCondition(geometry_context, 2, 3)
     
    lb1 = LineBundle(2*H, geometry_context=geometry_context)
    lb2 = LineBundle(3*H, geometry_context=geometry_context)

    ## From previous test we know that the central charge of lb1 is 17

    assert sigma1.phase(lb1) == 0.0

    ## From previous test we know that the central charge of lb2 is 15 + 12i

    assert sigma1.phase(lb2) == 0.21477671252272273

    gco = GradedCoproductObject([lb1], shift_vector=[4])
    print(sigma1.phase(gco))



def test_StabilityCondition_get_HN_factors():
    H, C, D = symbols("H C D")
    basis = [H]
    tensor_data = {
        (H, H): 4,
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)

    sigma1 = StabilityCondition(geometry_context, 2, 2)
    
    lb1 = LineBundle(H, geometry_context=geometry_context)
    lb2 = LineBundle(2*H, geometry_context=geometry_context)
    lb3 = LineBundle(3*H, geometry_context=geometry_context)
    lb4 = LineBundle(4*H, geometry_context=geometry_context)
    lb5 = LineBundle(5*H, geometry_context=geometry_context)

    sph=SphericalTwistComposition([lb3, lb1])


    print("\n\n")
    print(sph.defining_triangle)
    print("\n\n")


    print(sigma1.get_HN_factors(sph))



def test_StabilityCondition_K3_line_bundle_destab():


    H, C, D = symbols("H C D")
    basis = [H, C, D]
    tensor_data = {
        (H, H): 4,
        (C,C):-2,
        (H,C):1,
        (D,D):-2,
        (H,D):2,
        (C,D):0
    }

    divisor_data = DivisorData(basis=basis, top_intersection_form=tensor_data)
    geometry_context = GeometryContext(catagory='K3', divisor_data=divisor_data, polarization=H)

    print("\n\n")

    Num_iter = 200

    import random
    import time
    # from tqdm import tqdm

    start_time = time.perf_counter()
    

    # for _ in tqdm(range(Num_iter), desc="Iterations"):
    #     r1 = random.uniform(-10, 10)
    #     r2 = random.uniform(-10, 10)
    #     r3 = random.uniform(-10, 10)
    #     r4 = random.uniform(1, 10)

        

    #     C_coeff = random.randint(-10, 10)
    #     D_coeff = random.randint(-10, 10)
    #     H_coeff = random.randint(-10, 10)

    #     print("\n\n: Random parameters: ", r1, r2, r3, r4)
    #     print("Divisor: ", C_coeff, "*C + ", D_coeff, "*D + ", H_coeff, "*H\n\n")

    r1 = 2.034692616282385
    r2= 4.732391830282198
    r3=-8.053282369672061
    r4=4.3042163105497355
    C_coeff = 9
    D_coeff = 9
    H_coeff = 10

    print(f"Random parameters: {r1}, {r2}, {r3}, {r4}\n\n Divisor: {C_coeff}*C + {D_coeff}*D + {H_coeff}*H\n\n")



    sigma = StabilityCondition(geometry_context, r1,r2,r3,r4)

    lb1 = LineBundle(C_coeff*C + D_coeff*D + H_coeff*H, geometry_context=geometry_context)

    destabilizers = sigma._K3_line_bundle_destabilizing_candidates(lb1, r_max=3)
    if destabilizers is not None:
        print(f"Stability parameters: {r1}, {r2}, {r3}, {r4}")
        print(f"Line bundle: {C_coeff}*C + {D_coeff}*D + {H_coeff}*H")
        print("Destabilizers:")
        for d in destabilizers:
            print(d)

    end_time = time.perf_counter()
    print(f"Execution time for {Num_iter} calls: {end_time - start_time} seconds")
    print("\n\n")




    