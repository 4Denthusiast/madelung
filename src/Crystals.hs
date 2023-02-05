{-# LANGUAGE BlockArguments #-}

module Crystals (
) where

import Vectors

import Debug.Trace
import Numeric.SpecFunctions (erfc)

cubicLattice :: Int -> Lattice
cubicLattice n = map (\i -> replicate i 0 ++ [1] ++ replicate (n-i-1) 0) [0..(n-1)]

bccLattice :: Int -> Lattice
bccLattice n = replicate n 0.5 : tail (cubicLattice n)

fccLattice :: Int -> Lattice
fccLattice n = (1:(-1):replicate (n-2) 0) : map (1:) (cubicLattice (n-1))

-- Includes a simplex of the full dimension. Equivalent to fcc in 3D.
simplexLattice :: Int -> Lattice
simplexLattice n = map (\i -> replicate i x ++ [x'] ++ replicate (n-i-1) x) [0..(n-1)]
    where x = (sqrt (fromIntegral (1+n)) - 1) / (sqrt 2 * fromIntegral n)
          x' = x + sqrt 0.5

data Crystal = Crystal Lattice [(Double, Vect)] deriving Show

nacl :: Int -> Crystal
nacl n = Crystal (fccLattice n) [(1,replicate n 0),(-1,1:replicate (n-1) 0)]

cscl :: Int -> Crystal
cscl n = Crystal (cubicLattice n) [(1,replicate n 0),(-1,replicate n 0.5)]

zincblende :: Int -> Crystal
zincblende n = Crystal l [(1,replicate n 0),(-1,replicate n y)]
    where l = simplexLattice n
          y = 1/(2*sum (head l)) --equidistant from 0 and all the basis vectors

jellium :: Lattice -> Crystal
jellium l = Crystal l [(1,scale 0 (head l))]

-- The energy of the unit cell in the total field (without particles self-interacting)
totalEnergy :: Crystal -> Double
totalEnergy c@(Crystal l as) = energyWithCutoff c (det l ** (1/fromIntegral (length (head l))))

energyWithCutoff :: Crystal -> Double -> Double
energyWithCutoff c a = shortRangeEnergy c a + longRangeEnergy c a

-- The energy of the unit cell in the convolution of the field with a^-d exp(-pi r^2/a^2), including self-interaction, calculated in reciprocal space.
longRangeEnergy :: Crystal -> Double -> Double
longRangeEnergy (Crystal l as) a = sumOver baseField \(v,s) -> s * ((sumOver as \(q,v') -> q*cos(2*pi*dot v v'))^2 + (sumOver as \(q,v') -> q*sin(2*pi*dot v v'))^2)
    where d :: Num n => n
          d = fromIntegral $ dimension l
          -- The field, in reciprocal space, of a lattice l of gaussians with a compensatory uniform background.
          baseField = map (\v -> (v,(d-2)*sphereArea d/(4*pi*pi)/(dot v v)*exp(-pi*a*a*dot v v)/abs(det l))) $ filter (not . isZero) $ generateLattice (3/a) $ dual l

-- The energy of the unit cell in the total field minus the long range energy.
shortRangeEnergy :: Crystal -> Double -> Double
shortRangeEnergy (Crystal l as) a = backgroundTerm + sumOver (generateLattice (4*a) l) \v -> sumOver as \(q0,v0) -> sumOver as \(q1,v1) -> q0*q1*potential (dimension l) (norm (add v (add v0 (scale (-1) v1))))
    where potential d 0 = - sphereArea d * a^^(2-d) / (2*pi) -- Self interaction is ignored but it's included in the long-range, so it still needs to be cancelled.
          potential 3 r = erfc(r*sqrt pi/a)/r
          potential 4 r = exp(-pi*r*r/(a*a)) / (r*r)
          backgroundTerm = -(sumOver as fst)^2 * case (dimension l) of {3 -> a^2; 4 -> pi*a^2} / abs (det l)

-- The surface area of a unit sphere in d-dimensional space i.e. S^(d-1)
sphereArea :: Int -> Double
sphereArea d = (2*pi)^(div d 2)*2^(mod d 2) / doubleFac (d-2)
    where doubleFac x | x < 2 = 1
                      | otherwise = fromIntegral x * doubleFac (x-2)

sumOver :: Num n => [a] -> (a -> n) -> n
sumOver xs f = sum $ map f xs

minSeparation :: Crystal -> Double
minSeparation (Crystal l as) = minimum $ map norm $ filter (not . isZero) $ add <$> aos <*> generateLattice ((2*) $ maximum $ map norm (aos ++ l)) l
    where aps = map snd as
          aos = (add . scale (-1)) <$> aps <*> aps

madelungsConstant :: Crystal -> Double
madelungsConstant c@(Crystal l _) = totalEnergy c * minSeparation c ^(dimension l-2)

benchmark :: Int -> [Double]
benchmark d = map madelungsConstant $ [nacl d, cscl d, zincblende d] ++ map jellium [cubicLattice d, bccLattice d, fccLattice d, simplexLattice d]
