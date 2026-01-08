#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from geopy.geocoders import Nominatim
import pandas as pd
import time
import folium
import matplotlib.pyplot as plt


cities_states = [
    "Honolulu, HI", "Anaheim, CA", "Henderson, NV", "Orlando, FL", "Lexington, KY", "Stockton, CA",
    "Riverside, CA", "Corpus Christi, TX", "Irvine, CA", "Cincinnati, OH", "Santa Ana, CA", "Newark, NJ",
    "St. Paul, MN", "Pittsburgh, PA", "Greensboro, NC", "Durham, NC", "Lincoln, NE", "Jersey City, NJ",
    "Plano, TX", "Anchorage, AK", "North Las Vegas, NV", "St. Louis, MO", "Madison, WI", "Chandler, AZ",
    "Gilbert, AZ", "Reno, NV", "Buffalo, NY", "Chula Vista, CA", "Fort Wayne, IN", "Lubbock, TX",
    "Toledo, OH", "St. Petersburg, FL", "Laredo, TX", "Irving, TX", "Chesapeake, VA", "Glendale, AZ",
    "Winston-Salem, NC", "Port St. Lucie, FL", "Scottsdale, AZ", "Garland, TX", "Boise City, ID",
    "Norfolk, VA", "Spokane, WA", "Richmond, VA", "Fremont, CA", "Huntsville, AL", "Frisco, TX",
    "Cape Coral, FL", "Santa Clarita, CA", "San Bernardino, CA", "Tacoma, WA", "Hialeah, FL", "Baton Rouge, LA",
    "Modesto, CA", "Fontana, CA", "McKinney, TX", "Moreno Valley, CA", "Des Moines, IA", "Fayetteville, NC",
    "Salt Lake City, UT", "Yonkers, NY", "Worcester, MA", "Rochester, NY", "Sioux Falls, SD", "Little Rock, AR",
    "Amarillo, TX", "Tallahassee, FL", "Grand Prairie, TX", "Columbus, GA", "Augusta-Richmond County, GA",
    "Peoria, AZ", "Oxnard, CA", "Knoxville, TN", "Overland Park, KS", "Birmingham, AL", "Grand Rapids, MI",
    "Vancouver, WA", "Montgomery, AL", "Huntington Beach, CA", "Providence, RI", "Brownsville, TX", "Tempe, AZ",
    "Akron, OH", "Glendale, CA", "Chattanooga, TN", "Fort Lauderdale, FL", "Newport News, VA", "Mobile, AL",
    "Ontario, CA", "Clarksville, TN", "Cary, NC", "Elk Grove, CA", "Shreveport, LA", "Eugene, OR", "Aurora, IL",
    "Salem, OR", "Santa Rosa, CA", "Rancho Cucamonga, CA", "Pembroke Pines, FL", "Fort Collins, CO",
    "Springfield, MO", "Oceanside, CA", "Garden Grove, CA", "Lancaster, CA", "Murfreesboro, TN", "Palmdale, CA",
    "Corona, CA", "Killeen, TX", "Salinas, CA", "Roseville, CA", "Denton, TX", "Surprise, AZ", "Macon-Bibb County, GA",
    "Paterson, NJ", "Lakewood, CO", "Hayward, CA", "Charleston, SC", "Alexandria, VA", "Hollywood, FL",
    "Springfield, MA", "Kansas City, KS", "Sunnyvale, CA", "Bellevue, WA", "Joliet, IL", "Naperville, IL",
    "Escondido, CA", "Bridgeport, CT", "Savannah, GA", "Olathe, KS", "Mesquite, TX", "Pasadena, TX",
    "McAllen, TX", "Rockford, IL", "Gainesville, FL", "Syracuse, NY", "Pomona, CA", "Visalia, CA",
    "Thornton, CO", "Waco, TX", "Jackson, MS", "Columbia, SC", "Fullerton, CA", "Torrance, CA", "Victorville, CA",
    "Midland, TX", "Orange, CA"
]   #https://www.biggestuscities.com/, >138k & <360k, 2024 population

geolocator = Nominatim(user_agent="city_geocoder", timeout=10)

# Geocode each city and store coordinates
coordinates = {}
for city_state in cities_states:
    location = geolocator.geocode(city_state)
    if location:
        coordinates[city_state] = (location.latitude, location.longitude)
    time.sleep(2)  

df_coordinates = pd.DataFrame.from_dict(coordinates, orient='index', columns=['Latitude', 'Longitude'])

df_coordinates.to_csv("city_coordinates.csv", index=True)


df = pd.read_csv("city_coordinates.csv")

# Create a map centered at an average location
map_cities = folium.Map(location=[39.8283, -98.5795], zoom_start=4)

for _, row in df.iterrows():
    folium.Marker([row['Latitude'], row['Longitude']], popup=row['Unnamed: 0']).add_to(map_cities)

map_cities.save("htl_us_cities_map.html")
print("Interactive map saved as 'htl_us_cities_map.html'.")


plt.title("Tentative HTL Plant Locations")
plt.show()

