# t1 = ['V', 'A', 'I', 'G', 'U', 'M', 'K', 'X', 'P', 'O', 'N', 'D', 'C', 'Y', 'Q', 'F', 'R', 'S', 'Z', 'L', 'E', 'W', 'T', 'B', 'H']
# t2 = ['H', 'F', 'M', 'W', 'S', 'R', 'L', 'A', 'G', 'I', 'N', 'K', 'D', 'V', 'E', 'X', 'C', 'Q', 'T', 'P', 'Y']


t1 = ['V', 'A', 'C', 'G', 'U']
t2 = ['V', 'A', 'I', 'G', 'U', 'K']

set1 = set(t1)
set2 = set(t2)

t3 = set1 - set2
print(t3, len(t3))

# {'O', 'U', 'Z', 'B'} 4

for i in range(len(t2)):
    print(f"'{t2[i]}': {i+1},")