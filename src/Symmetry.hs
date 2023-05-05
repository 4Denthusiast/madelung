{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Symmetry (
) where

import Data.Complex
import Data.Function
import Data.List
import Data.Ord
import Data.Ratio
import Debug.Trace
import qualified Numeric.LinearAlgebra as V
import Util (nTimes)

import Vectors

-- For any symmetry class, there is some basis for the lattice in which every basis vector is sent to one of these points by each of the symmetries in the point group.
nearZero :: Int -> [[Int]]
nearZero d = (\s v -> map (*s) v) <$> [-1,1] <*> filter (not . all (==0)) (sequence $ replicate d [0,1])

instance Fractional Int where
    a / b = if mod a b == 0 then div a b else error (show a ++ " is not divisible by " ++ show b)
    fromRational r = fromIntegral (numerator r) / fromIntegral (denominator r)

-- Enumerate all the possible elements of all the point groups, expressed in a basis of lattice vectors (along with some spurious extras involving shearing and dilation).
allPointSymmetries :: Int -> [[[Int]]]
allPointSymmetries d = nTimes d (concatMap extend) [[]]
    where extend s = filter (not . singular) ((:s) <$> nearZero d)
singular [] = False
singular s = any (all (==0)) s || singular (reduce (sortBy (on divisibility head) s))
divisibility x y | x == y = EQ
               | x == 0 = GT
               | y == 0 = LT
               | otherwise = compare (abs x) (abs y)
reduce ((0:p):s) = p:map tail s
reduce ((ph:p):s) = map (\(vh:v) -> zipWith (+) v $ map (*(-vh/ph)) p) s

type IntMat = V.Matrix V.I

toIntMat :: [[Int]] -> IntMat
toIntMat = V.tr . V.fromLists . map (map fromIntegral) -- transpose because I want it as columns. I'm not sure whether this is actually necessary.

isCoprime :: Int -> Int -> Bool
isCoprime 0 n = False
isCoprime 1 n = True
isCoprime n m | n < 0 = isCoprime (-n) m
              | n > m = isCoprime m n
              | otherwise = isCoprime (mod m n) n

--All the matrices X such that M^n=X and X^n'=M for some n and n', only including those that are possible outputs from allPointSymmetries
cyclicEquivalents :: IntMat -> Maybe [IntMat]
cyclicEquivalents = fmap selectCoprimes . cycle
    where cycle :: IntMat -> Maybe [IntMat]
          cycle m = cycle' (V.ident (V.rows m)) m [] m
          cycle' i m mms mm | mm == i = Just (mm:mms)
                            | not $ generateable mm = Nothing
                            | otherwise = cycle' i m (mm:mms) (mm <> m)
          selectCoprimes l = map snd $ filter (isCoprime (length l) . fst) $ zip [0..] l

generateable :: IntMat -> Bool
generateable = all (\v -> zeroone v || zeroone (-v)) . V.toColumns
zeroone = all (\n -> n == 0 || n == 1) . V.toList

-- For each cyclic subgroup of a possible point group, a single arbitrary representative.
basicPointSymmetries :: Int -> [IntMat]
basicPointSymmetries = filter bestInCycle . map toIntMat . allPointSymmetries
    where bestInCycle m = case cyclicEquivalents m of
            Nothing -> False
            Just ms -> all (>=m) ms

-- Checks whether the matrix involves shearing, and is therefore not a valid element of the special linear group for any metric. If the matrix is defective, eig returns nearly parallel eigenvectors in an attempt to represent it. There are probably large (either in magnitude or dimension) matrices where this value of epsilon is insufficient, but it seems to be working fine for this use.
defective :: IntMat -> Bool
defective = (<epsilon) . magnitude . V.det . snd . V.eig . (V.fromInt :: IntMat -> V.Matrix V.C)
    where epsilon = 1e-3

-- All of the signed permutation matrices of the given order, with their inverses.
signedPermutationMatrices :: Int -> [(IntMat, IntMat)]
signedPermutationMatrices d = combine <$> permutationMatrices <*> signedMatrices
    where combine p s = (p <> s, s <> V.tr p)
          permutationMatrices = map V.fromLists $ permutations $ map (\i -> replicate i 0 ++ [1] ++ replicate (d-i-1) 0) [0..(d-1)]
          signedMatrices = map ((d V.>< d) . intercalate (replicate d 0) . map (:[])) $ sequence $ replicate d [-1,1]

-- One representative of each conjugacy class under the subgroup of signed permutation matrices (I don't have a principled reason to choose that over the conjugacy classes in the special linear group as a whole, that just seemed like it would take too long)
representativePointSymmetries :: Int -> [IntMat]
representativePointSymmetries d = foldr (filter . betterThanConjugate) (basicPointSymmetries d) (signedPermutationMatrices d)
    where betterThanConjugate (p,p') m = let c = p' <> m <> p in not (generateable c) || m <= c

-- The first argument is the generators of the group. The second is all the elements of the group. The third represents the restrictions on the lattice shape in some way I haven't entirely worked out yet.
data PointGroup = PG [IntMat] [IntMat] (V.Matrix Double)

-- May include some groups that are conjugate to each other.
allPointGroups :: Int -> [PointGroup]
allPointGroups d = removeExactDuplicates $ concatMap expand $ initialGroups
    where initialGroups = map (\s -> (PG [s] (addGroupElement [V.ident d] s) (normRestrictions s), map (\m -> (m,normRestrictions m)) (basicPointSymmetries d))) (representativePointSymmetries d)
          expand = expandTrimmed . trimRemainingSymmetries
          trimRemainingSymmetries (PG gens g rs, ss) = (PG gens g rs, filter (\(s,sn) -> head gens <= s && not (elem s g) && normValid d (V.orth (rs V.||| sn))) ss)
          expandTrimmed (g, ss) = g : concatMap (expandWith g ss) ss
          expandWith (PG gens g rs) ss (s,sn) = expand $ (PG (s:gens) (addGroupElement g s) (V.orth (rs V.||| sn)), ss)
          addGroupElement g e = addGroupElement' (e:g) g []
          addGroupElement' :: [IntMat] -> [IntMat] -> [IntMat] -> [IntMat]
          addGroupElement' gens new old = if new == [] then old else addGroupElement' gens (filter (flip notElem (new ++ old)) ((<>) <$> new <*> gens)) (new ++ old)
          removeExactDuplicates = map (snd . head) . groupBy (\(x,_) (y,_) -> x==y) . sortOn fst . map (\pg@(PG _ g _) -> (sort g, pg))

-- Considering the norm as a d(d+1)/2 dimensional vector (e1.e1, e2.e2, ..., ed.ed, e1.e2, e2.e3, ..., e1.ed) the transform m preserving the norm is equivalent to every column of (normRestrictions m) being perpendicular to it. Ideally I'd keep it as ints to avoid rounding errors given I'll need equality comparisons using the result, but the algorithm to compute the range of an int matrix is more complicated and not in the library.
normRestrictions :: IntMat -> V.Matrix Double
normRestrictions m = V.orth $ V.fromLists elements
    where elements = map (\i -> map (element i) (triangleIndices (V.rows m))) (triangleIndices (V.rows m))
          element (qx,qy) (ix,iy) = fromIntegral $ (if (ix,iy)==(qx,qy) then -1 else 0) + V.atIndex m (qx,ix) * V.atIndex m (qy,iy) + (if qx /= qy then V.atIndex m (qy,ix) * V.atIndex m (qx,iy) else 0)

triangleIndices d = [(x,y) | x <- [0..d-1], y <- [0..x]]

-- Is there a positive definite norm which satisfies these restrictions? I don't have a proof that the method I'm using is correct but it seems plausible. It searches iteratively for a matrix which is the outer product of a vector with itself and is in the span of the argument by repeatedly approximating the matrix with a single square by taking the square of the eigenvector corresponding to the greatest eigenvalue then projecting that onto the span of the argument. Because of the iteration this is terribly inefficient.
normValid :: Int -> V.Matrix Double -> Bool
normValid d rs = trace "------" $ V.cols rs == 0 || getsStuck (-1) (head $ V.toColumns $ traceShowId rs)
    where epsilon = 1e-2
          getsStuck pl v = traceShow v $ let
                  v' = traceShowId $ project $ traceShowId $ V.normalize $ traceShowId $ flattenSymMatrix $ traceShowId $ positiveEigenvectors $ traceShowId $ toSymmMatrix v
                  l = traceShowId $ V.norm_2 v'
                  v'' = V.scale (1/l) v'
              in (l >= 0) && not (l > 1-epsilon) && (abs (l-pl) < epsilon*epsilon || getsStuck l v'')
          toSymmMatrix v = V.trustSym $ (d V.>< d) [v V.! reverseIndex (max x y) (min x y) | x <- [0..d-1], y <- [0..d-1]]
          reverseIndex x y = div (x*(x+1)) 2 + y
          positiveEigenvectors = (\(xs,vs) -> vs <> V.diag (V.cmap (max 0) xs) <> V.tr vs) . V.eigSH
          --sqrtDiag = V.cmap (sqrt . abs) . V.takeDiag
          --outerSquare = (\x -> x <> V.tr x)
          flattenSymMatrix :: V.Matrix Double -> V.Vector Double
          flattenSymMatrix m = V.fromList $ map (V.atIndex m) (triangleIndices d)
          project = (rs V.#>) . (V.<# rs)

instance Ord (V.Matrix V.I) where
    compare m1 m2 = compare (V.toLists m1) (V.toLists m2)

--For debugging purposes.
order :: IntMat -> Int
order m = order' m  1
    where order' mm n | mm == i                             = n
                      | not $ null $ V.find ((>1) . abs) mm = 0
                      | otherwise                           = order' (mm <> m) (n+1)
          i = V.ident (V.rows m)

subCycleSignature :: IntMat -> [Int]
subCycleSignature = reverse . sort . map cOrder . V.toList . fst . V.eig . (V.fromInt :: IntMat -> V.Matrix V.C)
    where epsilon = 1e-8
          cOrder x = cOrder' x x 1
          cOrder' x xx n | magnitude (xx - 1) < 1e-3 = n
                         | magnitude xx > 2 || magnitude xx < 0.5 = 0
                         | otherwise = cOrder' x (x*xx) (n+1)

counts :: forall a. Ord a => [a] -> [(a,Int)]
counts = sort . foldr insert []
    where insert :: a -> [(a,Int)] -> [(a,Int)]
          insert x cs = case lookup x cs of
            Nothing -> (x,1):cs
            Just n -> (x,n+1):filter (/=(x,n)) cs
