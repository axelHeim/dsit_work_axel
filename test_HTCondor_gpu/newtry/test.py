import sys

print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

liste = sys.argv
if len(liste) == 2:
    print(liste[1])
    print(type(liste[1]))