import matlab.engine

eng = matlab.engine.start_matlab()

eng.addpath(r'/home/camarotti/Simulations/CoSimulation /CoSimIO Tutorials')

result = eng.my_function(5)
print(result)

# Check MATLAB version
print(eng.version())

a = eng.sqrt(4.0)
print(a)

tf = eng.isprime(37)
print(tf)