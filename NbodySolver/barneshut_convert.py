import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Converts Princeton file format to csv.')
    parser.add_argument("princetonfile")
    args = parser.parse_args()

    meters_per_au = 149597871000
    kg_per_earthmass = 5.972e24
    seconds_per_day = 60*60*24

    with open(args.princetonfile, 'r') as f:
        contents = f.read()
    lines = contents.splitlines()
    
    print int(lines[0]) # nbodies
    print 0.0
    for line in lines[2:]:
        L = line.split()
        x = float(L[0]) / meters_per_au
        y = float(L[1]) / meters_per_au
        vx = float(L[2]) * seconds_per_day / meters_per_au
        vy = float(L[3]) * seconds_per_day / meters_per_au
        m = float(L[4]) / kg_per_earthmass
        
        print x,y,vx,vy,m

if __name__ == '__main__':
    main()
