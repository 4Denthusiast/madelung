{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module Symmetry (
) where

import Control.Monad.Trans.State
import Data.Complex
import Data.Function
import Data.List
import Data.Maybe
import Data.Ord
import Data.Ratio
import qualified Data.Set as S
import Data.Set (Set, member, notMember, singleton, empty, size)
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
data PointGroup = PG [IntMat] (Set IntMat) (V.Matrix Double) deriving (Show)

-- May include some groups that are conjugate to each other.
allPointGroups :: Int -> [PointGroup]
allPointGroups d = snd $ execState (mapM_ expand initialGroups) (empty,[])
    where initialGroups = zipWith (\n s -> (n,PG [s] (addGroupElement (singleton $ V.ident d) s) (normRestrictions s), map (\m -> (m,normRestrictions m)) (basicPointSymmetries d))) [0..] (representativePointSymmetries d)
          expand (n,pg@(PG gens g rs), ss) = do
              (gs,gl) <- get
              if member g gs then
                  trace (show n ++ map (const ' ') gens ++ show (size g, size gs)) $ return ()
              else trace (show n ++ map (const '-') gens ++ show (size g, length ss)) $ do
                  let ss' = trimRemainingSymmetries gens g rs ss
                  put (S.insert g gs, pg : gl)
                  mapM_ (expandWith n pg ss') ss'
          trimRemainingSymmetries gens g rs ss = filter (\(s,sn) -> head gens < s && notMember s g && normValid d (orth (rs V.||| sn)) && all (<=s) (filter (flip member $ S.fromList $ map fst ss) $ (s<>) <$> S.toList g)) ss
          expandWith n (PG gens g rs) ss (s,sn) = expand (n, PG (s:gens) (addGroupElement g s) (orth (rs V.||| sn)), ss)
          addGroupElement g e = addGroupElement' (e:S.toList g) g empty
          addGroupElement' :: [IntMat] -> Set IntMat -> Set IntMat -> Set IntMat
          addGroupElement' gens new old = if S.null new then old else addGroupElement' gens (S.difference (S.fromList $ (<>) <$> S.toList new <*> gens) (S.union new old)) (S.union new old)
          orth m = if V.cols m == 0 then m else V.orth m

-- Considering the norm as a d^2 dimensional vector, the transform m preserving the norm is equivalent to every column of (normRestrictions m) being perpendicular to it. Ideally I'd keep it as ints to avoid rounding errors given I'll need equality comparisons using the result, but the algorithm to compute the range of an int matrix is more complicated and not in the library.
normRestrictions :: IntMat -> V.Matrix Double
normRestrictions m = V.orth $ V.fromLists elements
    where elements = map (\i -> map (element i) (filter (uncurry (>=)) indices)) indices
          element (qx,qy) (ix,iy) = (halfElement (qx,qy) (ix,iy) + halfElement (qy,qx) (ix,iy)) / 2
          halfElement (qx,qy) (ix,iy) = fromIntegral $ (if (ix,iy)==(qx,qy) then -1 else 0) + V.atIndex m (qx,ix) * V.atIndex m (qy,iy)
          indices = [(x,y) | y <- [0..V.rows m - 1], x <- [0..V.rows m - 1]]

-- Is there a positive definite norm which satisfies these restrictions? I don't have a proof that the method I'm using is correct but it seems plausible. It searches iteratively for a matrix which is the outer product of a vector with itself and is in the span of the argument by repeatedly approximating the matrix with a single square by taking the square of the eigenvector corresponding to the greatest eigenvalue then projecting that onto the span of the argument. Because of the iteration this is terribly inefficient.
normValid :: Int -> V.Matrix Double -> Bool
normValid d rs = checkInterleaved (findNorm $ V.flatten $ V.ident d) (findObstruction $ head $ V.toColumns $ traceShowId rs)
    where {-checkInterleaved (True:_) _ = True
          checkInterleaved (False:_) (True:_) = False
          checkInterleaved (False:x) (False:y) = checkInterleaved x y-}
          checkInterleaved xs _ = any id $ take 4 xs
          
          findObstruction v = traceShow prettyRS $ let
                  v' = traceShowId $ V.flatten $ traceShowId $ positiveEigenvectors $ traceShowId $ toSymmMatrix v
                  e = traceShowId $ V.norm_2 (v' - v)
                  d = traceShowId $ V.norm_2 $ traceShowId $ project (v' - v)
                  v'' = traceShowId $ V.normalize $ traceShowId $ project $ traceShowId $ v + V.scale ((e/d)^2) (v' - v)
              in seq v'' $ (e < 1e-6) : findObstruction v''
          prettyRS = map (V.reshape d) $ V.toColumns $ V.cmap (\x -> if abs x < 1e-14 then 0 else x) $ prettify rs
          toSymmMatrix = V.trustSym . V.reshape d
          positiveEigenvectors = (\v -> v <> V.tr v) . (V.?? (V.All, V.Take 1)) . snd . V.eigSH
          project = (rs V.#>) . (V.<# rs)
          
          findNorm v = let
                  (xs,vs) = V.eigSH $ toSymmMatrix $ v - project v
                  v' = V.normalize $ V.flatten $ vs <> V.diag (V.cmap (max (0.5/fromIntegral d)) xs) <> V.tr vs
              in (all (> 0.2/fromIntegral d) $ V.toList xs) : findNorm v'

-- Note! This does not check that the matrices are the same size. Doing that makes it run noticably slower.
instance Ord (V.Matrix V.I) where
    compare m1 m2 = compare (V.flatten m1) (V.flatten m2)

pointGroups :: Int -> [PointGroup]
pointGroups d = concatMap (nubBy equivalent . map snd) $ groupBy (on (==) fst) $ sortOn fst $ map (\p@(PG _ g _) -> ((groupSignature g, messiness g), p)) $ allPointGroups d
    where groupSignature = counts . map subCycleSignature . S.toList
          messiness = sum . map (V.sumElements . V.cmap abs) . S.toList
          -- equivalent actually checks the first is conjugate to a subgroup of the second, as it is assumed the groups are the same size.
          equivalent (PG gens _ _) (PG _ g _) = any (\(s',s) -> all (\gen -> member (s' <> gen <> s) g) gens) symmetriesAndInverses
          symmetriesAndInverses = mapMaybe (\s -> (, s) <$> inv s) $ map toIntMat $ allPointSymmetries d
          inv s = inv' s s
          inv' ss s
              | s <> ss == V.ident d  = Just ss
              | not (generateable ss) = Nothing
              | otherwise             = inv' (ss <> s) s

--For debugging purposes. The criterion (>10) is not principled.
order :: IntMat -> Int
order m = order' m  1
    where order' mm n | mm == i                             = n
                      | not $ null $ V.find ((>10) . abs) mm = 0
                      | otherwise                           = order' (mm <> m) (n+1)
          i = V.ident (V.rows m)

--The group element order, except when the orbit goes outside the range of allowed matrices it counts as 0 even if it would come back again.
strictOrder :: IntMat -> Int
strictOrder m = order' m 1
    where order' mm n | mm == i               = n
                      | not (generateable mm) = 0
                      | otherwise             = order' (mm <> m) (n+1)
          i = V.ident (V.rows m)

subCycleSignature :: IntMat -> [Int]
subCycleSignature = reverse . sort . map cOrder . V.toList . fst . V.eig . (V.fromInt :: IntMat -> V.Matrix V.C)
    where epsilon = 1e-8
          cOrder x = cOrder' x x 1
          cOrder' x xx n | magnitude (xx - 1) < 1e-3 = n
                         | magnitude xx > 2 || magnitude xx < 0.5 = 0
                         | otherwise = cOrder' x (x*xx) (n+1)

groupSignature :: PointGroup -> [([Int],Int)]
groupSignature (PG _ g _) = counts $ map subCycleSignature $ S.toList g

counts :: forall a. Ord a => [a] -> [(a,Int)]
counts = sort . foldr insert []
    where insert :: a -> [(a,Int)] -> [(a,Int)]
          insert x cs = case lookup x cs of
            Nothing -> (x,1):cs
            Just n -> (x,n+1):filter (/=(x,n)) cs

-- A sparser matrix with the same range.
prettify :: V.Matrix Double -> V.Matrix Double
prettify = V.fromColumns . map (V.cmap roundClose . rescale) . prettify' . V.toColumns
    where prettify' = prettify'' 0 0 . map (V.cmap (\x -> if abs x > epsilon then x else 0))
          prettify'' x y vs
              | x >= length vs = vs
              | y >= length vs = prettify'' (x+1) 0 vs
              | x == y = prettify'' x (y+1) vs
              | subset (vs !! x) (vs !! y) = traceShow (x,y) $ prettify' $ sparsen (vs !! x) (vs !! y) : (take y vs ++ drop (y+1) vs)
              | otherwise = prettify'' x (y+1) vs
          subset x y = all id $ zipWith (<=) (map (/=0) $ V.toList x) (map (/=0) $ V.toList y)
          sparsen :: V.Vector Double -> V.Vector Double -> V.Vector Double
          sparsen x y = let pivot = head $ traceShowId $ (\ps -> filter ((>maximum (map abs ps)*0.1).abs) ps) $ traceShowId $ filter ((>0).abs) $ traceShowId $ zipWith (/) (V.toList x) (V.toList y) in if pivot /= pivot then error (show (x,y)) else y - V.scale (1/pivot) x
          rescale :: V.Vector Double -> V.Vector Double
          rescale v = V.scale (1/maximumBy (comparing abs) (V.toList v)) v
          roundClose x = if abs (fromIntegral (round (x*100)) - x*100) > epsilon*100 then x else fromIntegral (round (x*100))/100
          epsilon = 1e-14
