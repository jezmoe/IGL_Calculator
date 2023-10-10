# Ideal Gas Law Calculator
# Last Updated 10/10/2023 @ 1:35 PM


from decimal import Decimal, getcontext

getcontext().prec = 28

def psi_to_atm(psi):
    return Decimal(psi) / Decimal("14.5038")

def bar_to_atm(bar):
    return Decimal(bar) / Decimal("1.01325")

def atm_to_psi(atm):
    return Decimal(atm) * Decimal("14.5038")

def atm_to_bar(atm):
    return Decimal(atm) * Decimal("1.01325")

def fahrenheit_to_kelvin(fahrenheit):
    celsius = (Decimal(fahrenheit) - 32) * 5 / 9
    kelvin = celsius + Decimal("273.15")
    return kelvin

def calculate_molar_mass(element_type):
    molar_masses = {
        "O2": Decimal("32"), 
        "H2": Decimal("2.016"), 
        "N2": Decimal("28.0134"), 
        "CO2": Decimal("44.0095"), 
        "HE": Decimal("4.002")
    }
    return molar_masses.get(element_type, Decimal("0"))

def solve_for_volume(P, n, T, pressure_unit):
    R = Decimal("0.08314473")
    if pressure_unit == "psi":
        P = psi_to_atm(P)
    elif pressure_unit == "bar":
        P = bar_to_atm(P)
    V_liters = (Decimal(n) * R * Decimal(T)) / P
    V_cubic_feet = V_liters / Decimal("28.3168")
    return V_liters, V_cubic_feet

def solve_for_temperature(P, n, V, pressure_unit):
    R = Decimal("0.08314473")
    if pressure_unit == "psi":
        P = psi_to_atm(P)
    elif pressure_unit == "bar":
        P = bar_to_atm(P)
    T_kelvin = (Decimal(P) * Decimal(V)) / (Decimal(n) * R)
    return T_kelvin

def solve_for_pressure(n, T, V):
    R = Decimal("0.08314473")
    P_atm = (Decimal(n) * R * Decimal(T)) / Decimal(V)
    return P_atm

def main():
    


    """Calculates the volume, pressure, or temperature of a gas using the Ideal Gas Law."""
    while True:
        print("\nIdeal Gas Law Calculator")
        print("Last updated: 10/9/2023 - Jan Rubido\n")
        print("Select what you want to solve for:")
        print("1. Pressure")
        print("2. Temperature")
        print("3. Volume")

        choice = input("Enter your choice (1, 2, or 3): ")

        if choice == "1":
            print("\nYou are solving for the Pressure.")
            print("                |                    ")
            print("                |                    ")
            print("                V                    \n")
            pressure_unit = input("Enter the desired output pressure unit (Atm, psi, bar): ").strip().lower()
            
            V = float(input("Enter the volume (in Liters): "))

            element_type = input("Enter the type of element (O2, H2, N2, CO2, or He): ").strip().upper()
            print(f"Element Type is: {element_type}")           
            molar_mass = calculate_molar_mass(element_type)
            print(f"Molar Mass for {element_type} is: {molar_mass}")
            mass_grams = round(float(input("Enter the mass (in grams): ")), 3)
            n = Decimal(mass_grams) / molar_mass

            
            R = 0.08314473
            temperature_unit = input("Enter the temperature unit (C, K, or F): ").strip().lower()
            if temperature_unit == "c":
                T_celsius = round(float(input("Enter the temperature (in °C): ")), 3)
                T_kelvin = T_celsius + 273.15
            elif temperature_unit == "k":
                T_kelvin = round(float(input("Enter the temperature (in K): ")), 3)
            elif temperature_unit == "f":
                T_fahrenheit = round(float(input("Enter the temperature (in °F): ")), 3)
                T_kelvin = fahrenheit_to_kelvin(T_fahrenheit)   

            P_atm = solve_for_pressure(n, T_kelvin, V)            
            if pressure_unit == "psi":
                P = atm_to_psi(P_atm)
            elif pressure_unit == "bar":
                P = atm_to_bar(P_atm)
            else:
                P = P_atm

            P = float(P)

            

            print(f"\n|   The pressure is approximately: {P:.3f} {pressure_unit}  |")
            print(f"| The number of moles for {element_type} is: {n:.3f} moles |")

            another = input("\nDo you want to solve for another? (yes/no): ").strip().lower()
            if another == "yes":
                continue  # This will send it back to the start of the outer loop
            else:
                break

        elif choice == "2":
            print("\nYou are solving for the Temperature.")
            print("                   |                    ")
            print("                   |                    ")
            print("                   V                    \n")
            pressure_unit = input("Enter the pressure unit (Atm, psi, bar): ").strip().lower()
            P = round(float(input(f"Enter the pressure (in {pressure_unit}): ")), 3)
            V = float(input("Enter the volume (in Liters): "))
            element_type = input("Enter the type of element (O2, H2, N2, CO2, or He): ").strip().upper()
            molar_mass = calculate_molar_mass(element_type)
            mass_grams = round(float(input("Enter the mass (in grams): ")), 3)
            n = (Decimal(mass_grams) * 1000 / molar_mass) / 1000


            
            T_kelvin = solve_for_temperature(P, n, V, pressure_unit)
            T_celsius = T_kelvin - Decimal('273.15')

            T_fahrenheit = (T_kelvin - Decimal('273.15')) * Decimal('9')/Decimal('5') + Decimal('32')

            print(f"\nThe number of moles for {element_type} is: {n:.3f} moles")
            print(f"\n| Approximate Temperatures: |")
            print(f"| Celsius: {T_celsius:.2f}°C          |")
            print(f"| Kelvin: {T_kelvin:.2f}K           |")
            print(f"| Fahrenheit: {T_fahrenheit:.2f}°F       |\n")
            another = input("\nDo you want to solve for another? (yes/no): ").strip().lower()
            if another == "yes":
                continue  # This will send it back to the start of the outer loop
            else:
                break

        elif choice == "3":
            print("\nYou are solving for the Volume.")
            print("                |                    ")
            print("                |                    ")
            print("                V                    \n")
            pressure_unit = input("Enter the pressure unit (Atm, psi, bar): ").strip().lower()
            P = round(float(input(f"Enter the pressure (in {pressure_unit}): ")), 3)

        if pressure_unit not in ["atm", "psi", "bar"]:
            print("Invalid pressure unit. Please enter 'Atm', 'psi', or 'bar'.")
            continue

        temperature_unit = input("Enter the temperature unit (C, K, or F): ").strip().lower()
        if temperature_unit == "c":
            T_celsius = round(float(input("Enter the temperature (in °C): ")), 3)
            T_kelvin = T_celsius + 273.15
        elif temperature_unit == "k":
            T_kelvin = round(float(input("Enter the temperature (in K): ")), 3)
        elif temperature_unit == "f":
            T_fahrenheit = round(float(input("Enter the temperature (in °F): ")), 3)
            T_kelvin = fahrenheit_to_kelvin(T_fahrenheit)

        element_type = input("Enter the type of element (O2, H2, N2, CO2, or He): ").strip().upper()
        molar_mass = calculate_molar_mass(element_type)
        mass_grams = round(float(input("Enter the mass (in grams): ")), 3)
        n = (Decimal(mass_grams) * 1000 / molar_mass) / 1000


        if choice == "3":
            V_liters, V_cubic_feet = solve_for_volume(P, n, T_kelvin, pressure_unit)
            print(f"\nThe volume is approximately {V_liters:.3f} Liters and {V_cubic_feet:.3f} cubic feet.")
            print(f"The number of moles for {element_type} is: {n:.3f} moles")

        another = input("\nDo you want to solve for another? (yes/no): ").strip().lower()
        if another == "yes":
            continue  # This will send it back to the start of the outer loop
        else:
            break

if __name__ == "__main__":
    main()
