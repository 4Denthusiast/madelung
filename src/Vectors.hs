module Vectors (
    Vect,
    scale,
    dot,
    norm,
    add,
    isZero,
    Lattice,
    dimension,
    simplifyLattice,
    generateLattice,
    det,
    dual
) where

import Data.List
import Util (nTimes)

import Debug.Trace

type Vect = [Double]

dot :: Vect -> Vect -> Double
dot xs1 xs2 = sum $ zipWith (*) xs1 xs2

isZero :: Vect -> Bool
isZero = all (==0)

norm :: Vect -> Double
norm v = sqrt $ dot v v

add :: Vect -> Vect -> Vect
add = zipWith (+)

-- The use of toRational is probably horribly inefficient, but this doesn't really need to be high performance.
scale :: (Real a) => a -> Vect -> Vect
scale s = map ((*) $ fromRational $ toRational s)

type Lattice = [Vect]

dimension :: Lattice -> Int
dimension = length . head

-- A reasonably small (not necessarily minimal) vector that's the sum of the given vector with a lattive vector.
simplifyByLattice :: Lattice -> Vect -> Vect
simplifyByLattice l v = foldr simplifyByVector v l
    where simplifyByVector lv v' | isZero lv = v'
                                 | otherwise = add v' (scale (negate (round (dot lv v' / dot lv lv))) lv)

-- A reasonably small (not necessarily minimal) equivalent representation of the same lattice.
simplifyLattice :: Lattice -> Lattice
simplifyLattice l = if l' == l then l' else simplifyLattice l'
    where l' = filter (not . isZero) $ nTimes (length l) simplifyFirst l
          simplifyFirst (v : vs) = vs ++ [simplifyByLattice vs v]

-- All the lattice points within the given radius (hopefully).
generateLattice :: Double -> Lattice -> Lattice
generateLattice r l = filter (\v -> dot v v <= r * r) $ map (foldr1 add) $ sequence $ map (\v -> let r' = ceiling (r / norm v) in map (flip scale v) [-r'..r']) $ simplifyLattice $ l

det :: [[Double]] -> Double
det [] = 1
det vs = detWithPivot $ foldr findPivot (last vs,[],1) (init vs)
    where findPivot v t@(p,vs,s) = if abs (head p) > abs (head v) then (p,v:vs,-s) else (v,p:vs,s)
          detWithPivot (ph:pt,vs,s) | ph==0 = 0
                                    | ph/=0 = s * ph * det (map (\(h:t) -> add t (scale (-h/ph) pt)) vs)

dual :: Lattice -> Lattice
dual l = map (zipWith (*) (cycle [1,-1])) $ zipWith scale (cycle [1/det l,-1/det l]) $ map (map det . semiminors . transpose) $ semiminors l
    where semiminors x = zipWith (++) (inits x) (tail $ tails x)
