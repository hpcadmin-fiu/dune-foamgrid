// Refine the grid uniformly
template <int dimworld>
void Dune::FoamGrid<dimworld>::globalRefine (int refCount)
{
  willCoarsen=false;

  if (maxLevel()+refCount<0)
    DUNE_THROW(GridError, "Grid has only " << maxLevel() << " levels. Cannot do "
                           << " globalRefine(" << refCount << ")");

  // The leafiterator is simply successively visiting all levels from fine to coarse
  // and just checking whether the isLeaf flag is set. Using it to identify
  // elements that need refinement would produce an endless loop, as the newly
  // added elements
  // would later be identified as elements that still need refinement.
  std::size_t oldLevels =entityImps_.size();
  typedef typename
    std::vector<tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                      std::list<FoamGridEntityImp<1,dimworld> >,
                      std::list<FoamGridEntityImp<2,dimworld> > > >::reverse_iterator
                LevelIterator;

  // Allocate space for the new levels. Thus the rend iterator will not get
  // invalid due to newly added levels.
  entityImps_.reserve(oldLevels+refCount);
  levelIndexSets_.reserve(oldLevels+refCount);
  LevelIterator level = entityImps_.rbegin();

  // Add tuples for the new levels.
  for (int i=0; i < refCount; ++i)
  {
    entityImps_.push_back(tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                                std::list<FoamGridEntityImp<1,dimworld> >,
                                std::list<FoamGridEntityImp<2,dimworld> > >());
    levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
  }

  if (refCount < 0)
  {
    if (globalRefined+refCount<0)
    {
      for (int i=refCount; i<0; ++i)
      {
        // Mark each leaf element for coarsening.
        typedef typename Traits::template Codim<0>::LeafIterator Iterator;
        for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
             elem != end; ++elem)
        {
          mark(-1,*elem);
        }

        // do the adaptation
        preAdapt();
        adapt();
        postAdapt();
      }

      globalRefined=0;
      return;
    }
    else
    {
      for (int i=refCount; i<0; ++i)
      {
        delete levelIndexSets_.back();
        levelIndexSets_.pop_back();
      }

      entityImps_.resize(oldLevels+refCount);

      // To be able to create the leaf level we need to set
      // the sons of the entities of maxlevel to null
      typename std::list<FoamGridEntityImp<0,dimworld> >::iterator vIt
        = Dune::get<0>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<0,dimworld> >::iterator vEndIt
        = Dune::get<0>(entityImps_[maxLevel()]).end();
      for (; vIt!=vEndIt; ++vIt)
        vIt->son_=nullptr;

      typename std::list<FoamGridEntityImp<1,dimworld> >::iterator edIt
        = Dune::get<1>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<1,dimworld> >::iterator edEndIt
        = Dune::get<1>(entityImps_[maxLevel()]).end();
      for (; edIt!=edEndIt; ++edIt)
      {
        edIt->sons_[0]=nullptr;
        edIt->sons_[1]=nullptr;
        edIt->nSons_=0;
      }

      typename std::list<FoamGridEntityImp<2,dimworld> >::iterator elIt
        = Dune::get<2>(entityImps_[maxLevel()]).begin();
      typename std::list<FoamGridEntityImp<2,dimworld> >::iterator elEndIt
        = Dune::get<2>(entityImps_[maxLevel()]).end();
      for (; elIt!=elEndIt; ++elIt)
      {
        elIt->sons_[0]=nullptr;
        elIt->sons_[1]=nullptr;
        elIt->sons_[2]=nullptr;
        elIt->sons_[3]=nullptr;
        elIt->nSons_=0;
      }
    }

  }
  else
  {
    // sanity check whether the above asssumption really holds.
    assert(&entityImps_[oldLevels-1]==&(*level));

    // Do the actual refinement.
    // We start with the finest level
    std::size_t levelIndex;
    for (levelIndex=oldLevels-1; level!=entityImps_.rend(); ++level, --levelIndex)
    {
      typedef typename std::list<FoamGridEntityImp<2,dimworld> >::iterator ElementIterator;
      bool foundLeaf=false;

      for (ElementIterator element=Dune::get<2>(*level).begin(); element != Dune::get<2>(*level).end(); ++element)
        if(element->isLeaf())
        {
          foundLeaf = true;
          dverb << "refining element " << &(*element) << std::endl;
          if (element->type().isTriangle())
            refineSimplexElement(*element, refCount);
// 	  if (element->type().isLine())
// 	    refineLineElement(*element,refCount);
          else
            DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
        }

      if (!foundLeaf)
        break;
    }

    for (levelIndex+=2;levelIndex!=entityImps_.size(); ++levelIndex)
      levelIndexSets_[levelIndex]->update(*this, levelIndex);
  }

  // Update the leaf indices
  leafIndexSet().update(*this);

  globalRefined=std::max(globalRefined+refCount,0);
  postAdapt();
}


//f Book-keeping routine to be called before adaptation
template <int dimworld>
bool Dune::FoamGrid<dimworld>::preAdapt()
{
  // Loop over all leaf entities and check whether they might be
  // coarsened. If there is one return true.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;
  int addLevels = 0;
  willCoarsen = false;

  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>(); elem != end; ++elem)
  {
    int mark=getMark(*elem);
    addLevels=std::max(addLevels, elem->level()+mark-maxLevel());

    if (mark<0)
    {
      // If this element is marked for coarsening, but another child
      // of this element's father is marked for refinement, then we
      // need to reset the marker to doNothing
      bool otherChildRefined=false;
      FoamGridEntityImp<2,dimworld>& father = *this->getRealImplementation(*elem).target_->father_;
      typedef typename array<FoamGridEntityImp<2,dimworld>*,4>::iterator ChildrenIter;
      for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
        otherChildRefined = otherChildRefined ||
                            (*child)->markState_==FoamGridEntityImp<2,dimworld>::REFINE;

      if (otherChildRefined)
      {
        for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
          if ((*child)->markState_==FoamGridEntityImp<2,dimworld>::COARSEN)
            (*child)->markState_=FoamGridEntityImp<2,dimworld>::DO_NOTHING;
      }
      else
        willCoarsen = willCoarsen || mark<0;
    }
  }

  if (addLevels)
  {
    entityImps_.push_back(tuple<std::list<FoamGridEntityImp<0,dimworld> >,
                                std::list<FoamGridEntityImp<1,dimworld> >,
                                std::list<FoamGridEntityImp<2,dimworld> > >());
    levelIndexSets_.push_back(new FoamGridLevelIndexSet<const FoamGrid >());
  }

  return willCoarsen;
}


// Triggers the grid refinement process
template <int dimworld>
bool Dune::FoamGrid<dimworld>::adapt()
{
  std::set<std::size_t> levelsChanged;
  bool haveRefined=false;

  // Loop over all leaf elements and refine/coarsen those that marked for it.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;
  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>();
       elem != end; ++elem)
  {
    int mark=getMark(*elem);
    if (mark>0)
    {
      // Refine Triangles
      if (elem->type().isTriangle())
      {
        refineSimplexElement(*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_), 1);
        levelsChanged.insert(elem->level()+1);
        haveRefined=true;
      }
      else
        DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
    }

    if (mark<0) // If simplex was allready treated by coarsenSimplex mark will be 0
    {
      // Coarsen triangles
      if (elem->type().isTriangle())
      {
        assert(elem->level());
        coarsenSimplexElement(*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_));
        levelsChanged.insert(elem->level());
      }
      else
        DUNE_THROW(NotImplemented, "Refinement only supported for triangles!");
    }
  }

  if (!willCoarsen)
    return haveRefined;

  typedef typename std::set<std::size_t>::const_reverse_iterator SIter;
  for (SIter level=levelsChanged.rbegin(); level!=levelsChanged.rend(); ++level)
  {
    // First delete the pointer
    erasePointersToEntities(Dune::get<2>(entityImps_[*level]));

    // Now delete the actual vertices.
    eraseVanishedEntities(Dune::get<0>(entityImps_[*level]));

    // erase vanished edges
    {
      typedef typename std::list<FoamGridEntityImp<1,dimworld> >::iterator EdgeIter;
      for (EdgeIter edge=Dune::get<1>(entityImps_[*level]).begin(),
           edgeEnd=Dune::get<1>(entityImps_[*level]).end();
           edge != Dune::get<1>(entityImps_[*level]).end();)
      {
        bool hasSameLevelElements=false;
        typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::const_iterator ElementIter;
        for (ElementIter elem=edge->elements_.begin();
             elem != edge->elements_.end(); ++elem)
          hasSameLevelElements = hasSameLevelElements || (*elem)->level()==*level;

        assert(edge->willVanish_!=hasSameLevelElements);

        if (!hasSameLevelElements)
        {
          assert(edge->willVanish_);
          // erase returns next position
          edge=Dune::get<1>(entityImps_[*level]).erase(edge);
        }
        else
          // increment
          ++edge;
      }
    }

    // And the elements
    eraseVanishedEntities(Dune::get<2>(entityImps_[*level]));

    if (Dune::get<0>(entityImps_[*level]).size())
    {
      assert(Dune::get<1>(entityImps_[*level]).size() &&
             Dune::get<2>(entityImps_[*level]).size());
      // Update the level indices.
      levelIndexSets_[*level]->update(*this, *level);
    }
    else
    {
      assert(!Dune::get<1>(entityImps_[*level]).size() &&
             !Dune::get<2>(entityImps_[*level]).size());
      if (static_cast<int>(*level)==maxLevel())
        entityImps_.pop_back();
    }
  }

  if (levelsChanged.size())
    // Update the leaf indices
    leafIndexSet().update(*this);
  globalRefined=0;

  return haveRefined;
}



// Clean up refinement markers
template <int dimworld>
void Dune::FoamGrid<dimworld>::postAdapt()
{
  willCoarsen=false;

  // Loop over all leaf entities and remove the isNew Marker.
  typedef typename Traits::template Codim<0>::LeafIterator Iterator;

  for (Iterator elem=this->leafbegin<0>(), end = this->leafend<0>(); elem != end; ++elem)
  {
    FoamGridEntityImp<2,dimworld>& element=*const_cast<FoamGridEntityImp<2,dimworld>*>(this->getRealImplementation(*elem).target_);
    element.isNew_=false;
    assert(!element.willVanish_);
    if (element.father_)
      element.father_->markState_=FoamGridEntityImp<2,dimworld>::DO_NOTHING;
  }
}


// Erases pointers in father elements to vanished entities of the element
template <int dimworld>
void Dune::FoamGrid<dimworld>::erasePointersToEntities(std::list<FoamGridEntityImp<2,dimworld> >& elements)
{
  typedef typename std::list<FoamGridEntityImp<2,dimworld> >::iterator EntityIterator;
  for(EntityIterator element=elements.begin();
      element != elements.end(); ++element)
  {
    if(element->willVanish_)
    {
      FoamGridEntityImp<2,dimworld>& father=*element->father_;

      for (unsigned int i=0; i<father.nSons_; i++)
        father.sons_[i]=nullptr;
      for (unsigned int i=0; i<father.corners(); i++)
        if (father.vertex_[i]->son_!=nullptr)
          if (father.vertex_[i]->son_->willVanish_)
            father.vertex_[i]->son_=nullptr;
      for (unsigned int i=0; i<father.corners(); i++)
        if (father.edges_[i]->sons_[0]!=nullptr
            && father.edges_[i]->sons_[0]->willVanish_)
          for (unsigned int j=0; j<2; j++)
          {
            assert(father.edges_[i]->sons_[j]!=nullptr);
            assert(father.edges_[i]->sons_[j]->willVanish_);
            father.edges_[i]->sons_[j]=nullptr;
          }
    }
  }
}

// Erase Entities from memory that vanished due to coarsening.
template <int dimworld>
template<int i>
void Dune::FoamGrid<dimworld>::eraseVanishedEntities(std::list<FoamGridEntityImp<i,dimworld> >& levelEntities)
{
  typedef typename std::list<FoamGridEntityImp<i,dimworld> >::iterator EntityIterator;
  for (EntityIterator entity=levelEntities.begin();
      entity != levelEntities.end();)
  {
    if(entity->willVanish_)
      entity=levelEntities.erase(entity);
    else
      ++entity;
  }
}


// Coarsen an Element
template <int dimworld>
void Dune::FoamGrid<dimworld>::coarsenSimplexElement(FoamGridEntityImp<2,dimworld>& element)
{
  // If we coarsen an element, this means that we erase all chidren of its father
  // to prevent inconsistencies.
  const FoamGridEntityImp<2,dimworld>& father = *(element.father_);

  // The edges that might be erased
  std::set<FoamGridEntityImp<1,dimworld>*> childEdges;

  // The vertices that might be erased
  std::set<FoamGridEntityImp<0,dimworld>*> childVertices;
  typedef typename array<FoamGridEntityImp<2,dimworld>*,4>::const_iterator ChildrenIter;
  for (ChildrenIter child=father.sons_.begin(); child != father.sons_.end(); ++child)
  {
    // Remember element for the actual deletion taking place later
    (*child)->markState_=FoamGridEntityImp<2,dimworld>::IS_COARSENED;
    (*child)->willVanish_=true;

    typedef typename array<FoamGridEntityImp<1,dimworld>*,3>::iterator EdgeIter;

    for (EdgeIter edge=(*child)->edges_.begin(); edge != (*child)->edges_.end(); ++edge)
    {
      // Remove references to elements that will be erased
      typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
                    ElementIter;
      for (ElementIter element = (*edge)->elements_.begin();
           element != (*edge)->elements_.end();
           ++element)
      {
        for (ChildrenIter child1=father.sons_.begin(); child1 !=
             father.sons_.end(); ++child1)
          if(*element==*child1)
          {
            // To prevent removing an element from a vector
            // we just overwrite the element with its father.
            // later we will check the edges and erase all
            // of them that do not have any elements on the same
            // level.
            *element=&father;
          }
      }
    }

    typedef FoamGridEntityImp<0,dimworld>** VertexIter;
    for (VertexIter vertex=(*child)->vertex_; vertex != (*child)->vertex_+(*child)->corners();
        ++vertex)
      childVertices.insert(*vertex);
    typedef typename array<FoamGridEntityImp<1,dimworld>*, 3>::iterator EdgeIter;
    for (EdgeIter edge=(*child)->edges_.begin(); edge!= (*child)->edges_.end();
        ++edge)
      childEdges.insert(*edge);
  }

  // Check whether those guys are really erased.
  // That is, we remove all vertices and edges that are part of a child element of one
  // of the neighbours of the father element
  typedef typename array<FoamGridEntityImp<1,dimworld>*,3>::const_iterator EdgeIter;

  for(EdgeIter edge=father.edges_.begin(); edge != father.edges_.end(); ++edge)
  {
    typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator NeighborIter;
    for (NeighborIter neighbor = (*edge)->elements_.begin();
         neighbor != (*edge)->elements_.end();
         ++neighbor)
    {
      assert((*neighbor)->level()<=father.level());
      if(*neighbor != &father && (*neighbor)->level()==father.level() &&
         !(*neighbor)->isLeaf())
      {
        // This is a real neighbor element on the same level
        // Check whether one of its children is marked for coarsening
        bool coarsened=false;
        for (ChildrenIter child=(*neighbor)->sons_.begin();
             child != (*neighbor)->sons_.end(); ++child)
          if((*child)->mightVanish())
          {
            coarsened=true;
            break;
          }

        if(!coarsened)
        {
          // Remove all entities that exist in the children of the neighbor
          // from the list
          for (ChildrenIter child=(*neighbor)->sons_.begin();
               child != (*neighbor)->sons_.end(); ++child)
          {
            typedef typename array<FoamGridEntityImp<1,dimworld>*,3>::iterator EdgeIter;
            for (EdgeIter edge=(*child)->edges_.begin(); edge != (*child)->edges_.end(); ++edge)
              childEdges.erase(*edge);
            typedef FoamGridEntityImp<0,dimworld>** VertexIter;
            for (VertexIter vertex=(*child)->vertex_;
                 vertex != (*child)->vertex_+(*child)->corners();
                 ++vertex)
              childVertices.erase(*vertex);
          }
        }
      }
    }
  }

  typedef typename std::set<FoamGridEntityImp<1,dimworld>*>::iterator SEdgeIter;
  for (SEdgeIter e=childEdges.begin(); e!=childEdges.end(); ++e)
    (*e)->willVanish_=true;
  typedef typename std::set<FoamGridEntityImp<0,dimworld>*>::iterator VertexIter;
  for (VertexIter v=childVertices.begin(); v!=childVertices.end(); ++v)
   (*v)->willVanish_=true;
}



// Refine one element
template <int dimworld>
void Dune::FoamGrid<dimworld>::refineSimplexElement(FoamGridEntityImp<2,dimworld>& element,
                                                    int refCount)
{
  if(refCount<0)
  {
    // We always remove all the children from the father.
    // Removing means:
    // 1. remove pointers in father element, edges and vertices
    // 2. Mark removed entities for deletion.
    DUNE_THROW(NotImplemented, "Coarsening not implemented yet");
    return;
  }


  // TODO: Currently we are assuming that only globalRefine is available
  // and therefore the next level is not yet present.
  // For real adaptivity some of the edges and vertices might already be present.
  // Therefore we need some kind of detection for this later on.

  unsigned int nextLevel=element.level()+1;


  std::cout << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size()
            << std::endl;
  std::cout << "Edges " << nextLevel << ": " <<Dune::get<1>(entityImps_[nextLevel]).size()
            << std::endl;
  std::cout << "Elements " << nextLevel << ": " << Dune::get<2>(entityImps_[nextLevel]).size()
            << std::endl;

  array<FoamGridEntityImp<0,dimworld>*, 6> nextLevelVertices;
  std::size_t vertexIndex=0;

  // create copies of the vertices of the element
  for(unsigned int c=0; c<element.corners(); ++c)
  {
    dverb<<"Processing vertex "<<element.vertex_[c]<<std::endl;
    if(element.vertex_[c]->son_==nullptr){
      // Not refined yet
      Dune::get<0>(entityImps_[nextLevel])
      .push_back(FoamGridEntityImp<0,dimworld>(nextLevel,
                                               element.vertex_[c]->pos_,
                                               element.vertex_[c]->id_));
      FoamGridEntityImp<0,dimworld>& newVertex =
      Dune::get<0>(entityImps_[nextLevel]).back();
      element.vertex_[c]->son_=&newVertex;
    }
    check_for_duplicates(nextLevelVertices, element.vertex_[c]->son_, vertexIndex);
    nextLevelVertices[vertexIndex++]=element.vertex_[c]->son_;
  }

  // create new vertices from edge-midpoints together with the new edges that
  // have a father
  typedef typename array<FoamGridEntityImp<1,dimworld>*, 3>::iterator EdgeIterator;

  array<FoamGridEntityImp<1,dimworld>*, 9> nextLevelEdges;
  std::size_t edgeIndex=0;
  const Dune::ReferenceElement<double,dimension>& refElement
    = Dune::ReferenceElements<double, dimension>::general(element.type());

  // I am just to dumb for a general edge to vertice mapping.
  // Therefore we just store it here
  array<std::pair<unsigned int,unsigned int>,3 > edgeVertexMapping;
  edgeVertexMapping[0]=std::make_pair(0,1);
  edgeVertexMapping[1]=std::make_pair(2,0);
  edgeVertexMapping[2]=std::make_pair(1,2);

  for(EdgeIterator edge=element.edges_.begin(); edge != element.edges_.end(); ++edge)
  {
    typedef FoamGridEntityImp<0,dimworld> FoamGridVertex;
    const FoamGridVertex* v0 = element.vertex_[refElement.subEntity(edgeIndex/2, 1, 0, 2)];
    const FoamGridVertex* v1 = element.vertex_[refElement.subEntity(edgeIndex/2, 1, 1, 2)];

    if(!(*edge)->nSons_)
    {
      // Not refined yet
      // Compute edge midpoint
      FieldVector<double, dimworld> midPoint;
      for(int dim=0; dim<dimworld;++dim)
        midPoint[dim]=((*edge)->vertex_[0]->pos_[dim]
        + (*edge)->vertex_[1]->pos_[dim]) /2.0;

      //create midpoint
      Dune::get<0>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<0,dimworld>(nextLevel, midPoint,
                                               freeIdCounter_[0]++));
      FoamGridEntityImp<0,dimworld>& midVertex =
        Dune::get<0>(entityImps_[nextLevel]).back();
      check_for_duplicates(nextLevelVertices, &midVertex, vertexIndex);
      nextLevelVertices[vertexIndex++]=&midVertex;

      // sanity check for DUNE numbering
      dvverb<<"edge "<<edgeIndex/2<<": "<<"("<<(*edge)->vertex_[0]->son_<<","<<(*edge)->vertex_[1]->son_<<") with father ("<<(*edge)->vertex_[0]<<","<<(*edge)->vertex_[1]<<")"<<std::endl;
      assert(v0->son_!=nullptr);
      assert(v1->son_!=nullptr);
      assert(v0->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].first] ||
        v0->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].second]);

      assert(v1->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].first] ||
        v1->son_ == nextLevelVertices[edgeVertexMapping[edgeIndex/2].second]);

      // create the edges and publish them in the father
      Dune::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[edgeVertexMapping[edgeIndex/2].first], &midVertex,
                                                 nextLevel, freeIdCounter_[1]++, *edge));
      (*edge)->sons_[0] = &Dune::get<1>(entityImps_[nextLevel]).back();
      ++((*edge)->nSons_);
      // Inherit the boundaryId
      (*edge)->sons_[0]->boundaryId_=(*edge)->boundaryId_;
      nextLevelEdges[edgeIndex++]= (*edge)->sons_[0];

      // Initialize the elements_ vector of the new edge
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*edge)->sons_[0]->elements_=(*edge)->elements_;

      assert((*edge)->vertex_[1]->son_!=nullptr);
      Dune::get<1>(entityImps_[nextLevel])
        .push_back(FoamGridEntityImp<1,dimworld>(&midVertex, nextLevelVertices[edgeVertexMapping[edgeIndex/2].second],
                                                 nextLevel, freeIdCounter_[1]++, *edge));
      (*edge)->sons_[1] = &Dune::get<1>(entityImps_[nextLevel]).back();
      ++((*edge)->nSons_);
      // Inherit the boundaryId
      (*edge)->sons_[1]->boundaryId_=(*edge)->boundaryId_;
      nextLevelEdges[edgeIndex++]= (*edge)->sons_[1];

      // Initialize the elements_ vector of the new edge
      // with that of the father. Later we will overwrite it
      // with the correct values.
      (*edge)->sons_[1]->elements_=(*edge)->elements_;

      (*edge)->nSons_=2;
    } else {
      // Edges do already exist. Just add its sons to nextLevelEdges
      // but make sure that the one containing vertex edgeIndex comes first
      if((*edge)->sons_[0]->vertex_[0]->id_ ==
        nextLevelVertices[edgeVertexMapping[edgeIndex/2].first]->id_ ||
        (*edge)->sons_[0]->vertex_[1]->id_ ==
        nextLevelVertices[edgeVertexMapping[edgeIndex/2].first]->id_)
      {
        nextLevelEdges[edgeIndex++]=(*edge)->sons_[0];
      nextLevelEdges[edgeIndex++]=(*edge)->sons_[1];
      }else{
        nextLevelEdges[edgeIndex++]=(*edge)->sons_[1];
        nextLevelEdges[edgeIndex++]=(*edge)->sons_[0];
      }
      if((*edge)->sons_[0]->vertex_[0]->id_!=(*edge)->vertex_[0]->id_ &&
         (*edge)->sons_[0]->vertex_[0]->id_!=(*edge)->vertex_[1]->id_)
      {
        //vertex 0 is the midpoint
        check_for_duplicates(nextLevelVertices, (*edge)->sons_[0]->vertex_[0], vertexIndex);
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0,dimworld>*>((*edge)->sons_[0]->vertex_[0]);
      } else {
        check_for_duplicates(nextLevelVertices, (*edge)->sons_[0]->vertex_[1], vertexIndex);
        nextLevelVertices[vertexIndex++]=const_cast<FoamGridEntityImp<0,dimworld>*>((*edge)->sons_[0]->vertex_[1]);
      }
    }
  }
  assert(edgeIndex==6);
  // Create the edges that lie within the father element
  // first the one that lies opposite to the vertex 0 in the father
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[4], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  // the one opposite to father vertex 1
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[3],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  // and the one opposite to father vertex 2
  Dune::get<1>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<1,dimworld>(nextLevelVertices[4],
                                             nextLevelVertices[5], nextLevel,
                                             freeIdCounter_[1]++));
  nextLevelEdges[edgeIndex++]=&Dune::get<1>(entityImps_[nextLevel]).back();

  assert(edgeIndex==nextLevelEdges.size());
  assert(vertexIndex==nextLevelVertices.size());

  array<FoamGridEntityImp<2,dimworld>*, 4> nextLevelElements;
  // create the new triangles that lie in the corners
  // First the one that contains vertex 0 of the father.
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));

  FoamGridEntityImp<2,dimworld>* newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->edges_[0]=nextLevelEdges[0];
  newElement->edges_[1]=nextLevelEdges[3];
  newElement->edges_[2]=nextLevelEdges[6];
  newElement->vertex_[0]=nextLevelVertices[0];
  newElement->vertex_[1]=nextLevelVertices[3];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=0;
  nextLevelElements[0]=newElement;
  element.sons_[0]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;

  // Next the one that contains vertex 1 of the father.
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->edges_[0]=nextLevelEdges[4];
  newElement->edges_[1]=nextLevelEdges[1];
  newElement->edges_[2]=nextLevelEdges[7];
  newElement->vertex_[0]=nextLevelVertices[1];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[3];
  newElement->refinementIndex_=1;
  nextLevelElements[1]=newElement;
  element.sons_[1]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // Last the one that contains vertex 2 of the father.
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->edges_[0]=nextLevelEdges[2];
  newElement->edges_[1]=nextLevelEdges[5];
  newElement->edges_[2]=nextLevelEdges[8];
  newElement->vertex_[0]=nextLevelVertices[2];
  newElement->vertex_[1]=nextLevelVertices[4];
  newElement->vertex_[2]=nextLevelVertices[5];
  newElement->refinementIndex_=2;
  nextLevelElements[2]=newElement;
  element.sons_[2]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // create the triangle in the center
  Dune::get<2>(entityImps_[nextLevel])
    .push_back(FoamGridEntityImp<2,dimworld>(nextLevel, freeIdCounter_[2]++));
  newElement = &(Dune::get<2>(entityImps_[nextLevel]).back());
  newElement->isNew_=true;
  newElement->father_=&element;
  newElement->edges_[0]=nextLevelEdges[7];
  newElement->edges_[1]=nextLevelEdges[6];
  newElement->edges_[2]=nextLevelEdges[8];
  newElement->vertex_[0]=nextLevelVertices[3];
  newElement->vertex_[1]=nextLevelVertices[5];
  newElement->vertex_[2]=nextLevelVertices[4];
  newElement->refinementIndex_=3;
  nextLevelElements[3]=newElement;
  element.sons_[3]=newElement;
  dvverb<<"Pushed element "<<newElement<<" refindex="<<newElement->refinementIndex_<<std::endl;


  // Now that all the triangle are created, we can update the elements attached
  // to the edges.
  // The new neighbors of the edges lying on edges of the father element.
  std::size_t neighbors[6] = {0, 1, 2, 0, 1, 2};
  dvverb<<" element "<<&element<<std::endl;
  for(std::size_t i=0; i<6; ++i){
    // Overwrite the father element by the newly created elements.
    dvverb<<" neighbour "<<i<<": ";
    overwriteFineLevelNeighbours(*nextLevelEdges[i], nextLevelElements[neighbors[i]],
                                 &element);
  }

  // Update the neighbours of the inner edges
  nextLevelEdges[6]->elements_.push_back(nextLevelElements[0]);
  nextLevelEdges[6]->elements_.push_back(nextLevelElements[3]);
  nextLevelEdges[7]->elements_.push_back(nextLevelElements[3]);
  nextLevelEdges[7]->elements_.push_back(nextLevelElements[1]);
  nextLevelEdges[8]->elements_.push_back(nextLevelElements[3]);
  nextLevelEdges[8]->elements_.push_back(nextLevelElements[2]);

#ifndef NDEBUG
  for(std::size_t vidx=0; vidx<4; ++vidx)
  {
    dvverb<<std::endl<<" Element "<<nextLevelElements[vidx]<<": ";
    for(EdgeIterator edge=nextLevelElements[vidx]->edges_.begin();
        edge != nextLevelElements[vidx]->edges_.end(); ++edge)
    {
      dvverb<<std::endl<<"   edge "<<*edge<<": ";
      typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator
      ElementIterator;
      bool selfFound=false;
      for(ElementIterator elem=(*edge)->elements_.begin();
          elem != (*edge)->elements_.end();
          ++elem)
      {
        dvverb << *elem<<" ";
        if(*elem == nextLevelElements[vidx])
          selfFound=true;
      }
      assert(selfFound);
    }
  }

#endif
  element.nSons_=4;

  if((refCount--)>1)
  {
    int i=0;
    typedef typename array<FoamGridEntityImp<2,dimworld>*, 4>::iterator ElementIterator;
    for(ElementIterator elem=nextLevelElements.begin();
        elem != nextLevelElements.end(); ++elem)
    {
      dvverb<<std::endl<<"Refinining "<<(*elem)<<" (son of"<<&element<<") refCount="<<refCount<<" child="<<i++<<std::endl;
      refineSimplexElement(**elem, refCount);
    }
  }

  std::cout << "end refineSimplex" << std::endl;
  std::cout << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size()
            << std::endl;
  std::cout << "Edges " << nextLevel << ": " << Dune::get<1>(entityImps_[nextLevel]).size()
            << std::endl;
  std::cout << "Elements " << nextLevel << ": " << Dune::get<2>(entityImps_[nextLevel]).size()
            << std::endl;
}


// Overwrites the neighbours of this and descendant edges
template <int dimworld>
void Dune::FoamGrid<dimworld>::overwriteFineLevelNeighbours(FoamGridEntityImp<1,dimworld>& edge,
                                                            FoamGridEntityImp<2,dimworld>* son,
                                                            FoamGridEntityImp<2,dimworld>* father)
{
  typedef typename std::vector<const FoamGridEntityImp<2,dimworld>*>::iterator ElementIterator;
#ifndef NDEBUG
  bool fatherFound=false;
#endif
  for(ElementIterator elem=edge.elements_.begin();
      elem != edge.elements_.end();
      ++elem)
  {
    dvverb << *elem<<" ";
    if (*elem == father)
    {
#ifndef NDEBUG
      fatherFound=true;
#endif
      *elem = son;
    }
  }

#ifndef NDEBUG
  dvverb<<std::endl;
  assert(fatherFound);
#endif

  for(std::size_t i=0; i<edge.nSons_; ++i)
    overwriteFineLevelNeighbours(*edge.sons_[i], son, father);
}


// Recompute the grid indices after the grid has changed
template <int dimworld>
void Dune::FoamGrid<dimworld>::setIndices()
{
  // //////////////////////////////////////////
  //   Create the index sets
  // //////////////////////////////////////////
  for (int i=levelIndexSets_.size(); i<=maxLevel(); i++) {
    FoamGridLevelIndexSet<const FoamGrid >* p
      = new FoamGridLevelIndexSet<const FoamGrid >();
    levelIndexSets_.push_back(p);
  }

  for (int i=0; i<=maxLevel(); i++)
    if (levelIndexSets_[i])
      levelIndexSets_[i]->update(*this, i);

  // Update the leaf indices
  leafGridView_.indexSet_.update(*this);
}

//template <int dimworld>
//void Dune::FoamGrid<dimworld>::refineLineElement(FoamGridEntityImp<1,dimworld>& element,
//						 int refCount)
//{
//  std::cout<<"refineLineElement "<<std::endl;
//     if(refCount<0)
//   {
//     // We always remove all the children from the father.
//     // Removing means:
//     // 1. remove pointers in father element, edges and vertices
//     // 2. Mark removed entities for deletion.
//     DUNE_THROW(NotImplemented, "Coarsening not implemented yet");
//     return;
//   }
//
//   unsigned int nextLevel=element.level()+1;
//
//     std::cout << "Vertices " << nextLevel << ": " << Dune::get<0>(entityImps_[nextLevel]).size()
//               << std::endl;
//     std::cout << "Elements " << nextLevel << ": " <<Dune::get<1>(entityImps_[nextLevel]).size()
//               << std::endl;
//
//   array<FoamGridEntityImp<0,dimworld>*, 3> nextLevelVertices;
//   std::size_t vertexIndex=0;
//
//     for(unsigned int c=0; c<element.vertex_.size(); ++c)
//   {
//     //ES: kopieren der Punkte warum brauche ich das ?
//     if(element.vertex_[c]->son_==nullptr){
//             // Not refined yet
//       Dune::get<0>(entityImps_[nextLevel])
//       .push_back(FoamGridEntityImp<0,dimworld>(nextLevel,
//                                                element.vertex_[c]->pos_,
//                                                element.vertex_[c]->id_));
//       FoamGridEntityImp<0,dimworld>& newVertex =
//       Dune::get<0>(entityImps_[nextLevel]).back();
//       element.vertex_[c]->son_=&newVertex;
//
//     }
//   }
//
//       if(!(element->nSons_)
//     {
//       // Not refined yet
//       // Compute edge midpoint
//       FieldVector<double, dimworld> midPoint;
//       for(int dim=0; dim<dimworld;++dim)
//         midPoint[dim]=(element->vertex_[0]->pos_[dim]
//         + element->vertex_[1]->pos_[dim]) /2.0;
//
//       //create midpoint
//       Dune::get<0>(entityImps_[nextLevel])
//         .push_back(FoamGridEntityImp<0,dimworld>(nextLevel, midPoint,
//                                                freeIdCounter_[0]++));
//       FoamGridEntityImp<0,dimworld>& midVertex =
//         Dune::get<0>(entityImps_[nextLevel]).back();
//       check_for_duplicates(nextLevelVertices, &midVertex, vertexIndex);
//       nextLevelVertices[vertexIndex++]=&midVertex;
//     }
//
//}