"""
    trun_cone_area(r,R)

Function that calculates the area and volume of a truncated cone, given both radii.
The truncated cone represents our coffee cup.
"""
function trun_cone_area(r,R,h)
    slant_height = sqrt((R - r)^2 + h^2)
    lateral = π * (R + r) * slant_height
    top = π * r^2
    bottom = π * R^2

    SA = lateral + top + bottom
    V = (1/3) * π * h * (r^2 + r * R + R^2)
    print("surface area-")
    display(SA)
    print("volume-")
    display(V)

    return (SA, V)
end

"""
    dairy_percentage(water, dairy)

Calculates the overall heat coefficient of a water/dairy mixture based
on the proportions of each.
"""
function dairy_percentage(water_percent, water_density, water_heat, dairy_percent, dairy_density, dairy_heat, volume, SA)

    water_mass = water_density / volume * water_percent
    dairy_mass = dairy_density / volume * dairy_percent

    total_mass = water_mass + dairy_mass

    total_heat_capacity = ((water_mass / total_mass) * water_heat) + ((dairy_mass / total_mass) * dairy_heat)

    print("total heat capacity-")
    display(total_heat_capacity)

    thermal_cond = total_heat_capacity / (SA*(140-68))
    print("coffee/milk thermal conductivity-")
    display(thermal_cond)

    return thermal_cond
end

"""
    heat_transfer_co()

Function that calculates a 3-layer heat transfer coefficient.
Units for space in this function are milimeters.

Accepts a coffe temperature and a room temperature as variables. (In F)
"""
function heat_transfer_co(coffee_temp, room_temp)
    (SA, V) = trun_cone_area(25, 45, 100) # surface area of our truncated cone coffee cup in mm
    thick = 3 # thickness of the polypropelene cup

    dt = (coffee_temp - room_temp)

    air = 0.020 # thermal conductivity of air in 5mm layer around cone
    air_thermal = (10) / (air * SA)

    ceramic = 1.4 # thermal conductivity of ceramic
    ceramic_thermal = (2*thick) / (ceramic * SA)

    poly = 0.1 # 0.1 - 0.22 thermal conductivity
    poly_thermal = (5) / (poly * SA)

    water = 0.606 # thermal conductivity of water
    water_thermal = (5) / (water * SA)

    coffee = dairy_percentage(0.95, 1, 4.186, 0.05, 1.010, 3.62, V, SA) # thermal conductivity of water
    coffee_thermal = (5) / (coffee * SA)

    thermal_total = air_thermal + poly_thermal + coffee_thermal

    Q = (dt) / thermal_total

    print("overall heat transfer coefficient-")
    display(Q)
    return Q

end


ans = heat_transfer_co(140, 68)
display(ans)


"""

The overall heat transfer coefficient in this kinda of simulation takes a little bit of effort to calculate.
Examining the problem at hand, the overall heat transfer coefficient is an example of a three layered wall
heat transfer calculation. Heat transfer in a layered situation is given by the equation:

q = T_hot - T_cold / R

R is the sum of the heat transfer coefficients for each layer. In our simulation, we have three different layers total.
The coffee/milk mixture, the ceramic or polypropelene of the cup, and the air around the cup itself.

Water has a heat capacity of 4.186 J/g °C. Whole milk has a heat capacity of 3.89 J/g °C, while cream has a heat capacity of 3.35 J/g °C.
For a rough estimation of heat capacity of half and half, we averaged the heat capacity between whole milk and cream.


https://www.journalofdairyscience.org/article/S0022-0302(16)30053-4/pdf
https://www.researchgate.net/profile/Charles-Marzzacco/publication/331895987_Heat_Flow_in_a_Coffee_Cup/links/5c9238a4299bf1116939c7c7/Heat-Flow-in-a-Coffee-Cup.pdf
https://www.engineeringtoolbox.com/overall-heat-transfer-coefficient-d_434.html


"""
