module Symmetry (
) where

import Data.Function
import Data.List
import Data.Ratio
import Numeric.LinearAlgebra as V
import Util (nTimes)

import Vectors

-- For any symmetry class, there is some basis for the lattice in which every basis vector is sent to one of these points by each of the symmetries in the point group.
nearZero :: Int -> [[Int]]
nearZero d = (\s v -> map (*s) v) <$> [-1,1] <*> filter (not . all (==0)) (sequence $ replicate d [0,1])

instance Fractional Int where
    a / b = if mod a b == 0 then div a b else error (show a ++ " is not divisible by " ++ show b)
    fromRational r = fromIntegral (numerator r) / fromIntegral (denominator r)

-- Enumerate all the possible elements of all the point groups, expressed in a basis of lattice vectors.
allPointSymmetries :: Int -> [[[Int]]]
allPointSymmetries d = nTimes d (concatMap extend) [[]]
    where extend s = filter (not . singular) ((:s) <$> nearZero d)
          singular [] = False
          singular s = any (all (==0)) s || singular (reduce (sortBy (on divisibility head) s))
          divisibility x y | x == y = EQ
                           | x == 0 = GT
                           | y == 0 = LT
                           | otherwise = compare (abs x) (abs y)
          reduce ((0:p):s) = map tail s
          reduce ((ph:p):s) = map (\(vh:v) -> zipWith (+) v $ map (*(-vh/ph)) p) s
