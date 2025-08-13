def Function():
    data = [
        {"Name": "Sven", "Room number": "1.14", "Favorite food": "Lyoner"},
        {"Name": "Jens", "Room number": "1.13", "Favorite food": "Drohschele"},
        {"Name": "Johanna", "Room number": "1.13", "Favorite food": "Grumbeerkieschelscher"}
    ]

    print("Name,Room number,Favorite food")
    for i in data:
        print(f"{i['Name']},{i['Room number']},{i['Favorite food']}")

Function()
