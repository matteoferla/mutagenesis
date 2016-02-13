__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

class test:
    def __init__(self,data):
        self.data=data
    def __getitem__(self,x):
        print(x)
        return self.data[x[0]][x[1]]


class outer:
    class inner:
        def __init__(self,oself,a):
            self.y=oself.n+" from "+a

        def __str__(self):
            return self.y

    def __init__(self):
        self.n="Bye"
        self.x=self.inner(self,"Bob")


if __name__ == "__main__":
    print(outer().x)