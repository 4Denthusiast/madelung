{-# LANGUAGE BlockArguments #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}

module Crystals (
) where

import Vectors

import Debug.Trace
import GHC.TypeLits
import qualified Numeric.LinearAlgebra.Static as V
import Numeric.SpecFunctions (erfc)

cubicLattice :: KnownNat n => Lattice n
cubicLattice = l --This should just be V.eye, but that function is very buggy.
    where l = V.matrix $ concatMap (\i -> replicate i 0 ++ [1] ++ replicate (n-i-1) 0) [0..(n-1)]
          n = fromIntegral $ natVal l

bccLattice :: Dimension n => Lattice n
bccLattice = (cubicLattice V.=== V.row 0) V.||| V.col 0.5

fccLattice :: Dimension n => Lattice n
fccLattice = l
    where l = (cubicLattice V.||| V.col (V.vector (-1:replicate (fromIntegral (natVal l)-2) 0))) V.=== V.row 1

-- Includes a simplex of the full dimension. Equivalent to fcc in 3D.
simplexLattice :: Dimension n => Lattice n
simplexLattice = l
    where l = V.konst (sqrt 0.5) * cubicLattice + V.konst x
          x = (sqrt (fromIntegral (1+n)) - 1) / (sqrt 2 * fromIntegral n)
          n = fromIntegral $ natVal l

data Crystal n = Crystal (Lattice n) [(Double, V.R n)] deriving Show

nacl :: Dimension n => Crystal n
nacl = Crystal fccLattice [(1,0),(-1,0 V.& 1)]

cscl :: Dimension n => Crystal n
cscl = Crystal cubicLattice [(1,0),(-1,0.5)]

zincblende :: Dimension n => Crystal n
zincblende = Crystal l [(1,0),(-1,V.konst y)]
    where l = simplexLattice
          y = 1/(2*V.dot 1 (V.uncol $ fst $ V.splitCols l)) --equidistant from 0 and all the basis vectors

jellium :: Lattice n -> Crystal n
jellium l = Crystal l [(1,0)]

-- The energy of the unit cell in the total field (without particles self-interacting)
totalEnergy :: Dimension n => Crystal n -> Double
totalEnergy c@(Crystal l as) = energyWithCutoff c (abs(V.det l) ** (1/fromIntegral (natVal l)))

energyWithCutoff :: Dimension n => Crystal n -> Double -> Double
energyWithCutoff c a = shortRangeEnergy c a + longRangeEnergy c a

-- The energy of the unit cell in the convolution of the field with a^-d exp(-pi r^2/a^2), including self-interaction, calculated in reciprocal space.
longRangeEnergy :: Dimension n => Crystal n -> Double -> Double
longRangeEnergy (Crystal l as) a = sumOver baseField \(v,s) -> s * ((sumOver as \(q,v') -> q*cos(2*pi*V.dot v v'))^2 + (sumOver as \(q,v') -> q*sin(2*pi*V.dot v v'))^2)
    where d :: Num n => n
          d = fromIntegral $ natVal l
          -- The field, in reciprocal space, of a lattice l of gaussians with a compensatory uniform background.
          baseField = map (\v -> (v,(d-2)*sphereArea d/(4*pi*pi)/(V.dot v v)*exp(-pi*a*a*V.dot v v)/abs(V.det l))) $ filter ((/= 0).V.norm_1) $ generateLattice (3/a) $ dual l

-- The energy of the unit cell in the total field minus the long range energy.
shortRangeEnergy :: Dimension n => Crystal n -> Double -> Double
shortRangeEnergy (Crystal l as) a = backgroundTerm + sumOver (generateLattice (4*a) l) \v -> sumOver as \(q0,v0) -> sumOver as \(q1,v1) -> q0*q1*potential (fromIntegral $ natVal l) (V.norm_2 (v + v0 - v1))
    where potential d 0 = - sphereArea d * a^^(2-d) / (2*pi) -- Self interaction is ignored but it's included in the long-range, so it still needs to be cancelled.
          potential 3 r = erfc(r*sqrt pi/a)/r
          potential 4 r = exp(-pi*r*r/(a*a)) / (r*r)
          backgroundTerm = -(sumOver as fst)^2 * case (natVal l) of {3 -> a^2; 4 -> pi*a^2} / abs (V.det l)

-- The surface area of a unit sphere in d-dimensional space i.e. S^(d-1)
sphereArea :: Int -> Double
sphereArea d = (2*pi)^(div d 2)*2^(mod d 2) / doubleFac (d-2)
    where doubleFac x | x < 2 = 1
                      | otherwise = fromIntegral x * doubleFac (x-2)

sumOver :: Num n => [a] -> (a -> n) -> n
sumOver xs f = sum $ map f xs

minSeparation :: Dimension n => Crystal n -> Double
minSeparation (Crystal l as) = minimum $ filter (/=0) $ map V.norm_2 $ (+) <$> aos <*> generateLattice ((2*) $ maximum $ V.norm_2 l : map V.norm_2 aos) l
    where aps = map snd as
          aos = (-) <$> aps <*> aps

madelungsConstant :: Dimension n => Crystal n -> Double
madelungsConstant c@(Crystal l _) = totalEnergy c * minSeparation c ^(fromIntegral (natVal l)-2)

benchmark :: forall proxy n. Dimension n => proxy n -> [Double]
benchmark d = map madelungsConstant $ ([nacl, cscl, zincblende] ++ map jellium [cubicLattice, bccLattice, fccLattice, simplexLattice] :: [Crystal n])
