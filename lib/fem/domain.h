/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_DOMAIN_H
#define MECHSYS_DOMAIN_H

// Std Lib
#include <iostream>
#include <sstream>  // for istringstream, ostringstream
#include <cstdarg>  // for va_list, va_start, va_end
#include <map>

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/models/model.h>
#include <mechsys/fem/element.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/draw.h>

namespace FEM
{

typedef std::map<Node*,double>          InclSupport_t; ///< Inclined support type. Maps Node ==> alpha
typedef std::map<Node*,Array<double> >  Res_t;         ///< Maps node to results at nodes
typedef std::map<int, Array<Element*> > Tag2Eles_t;    ///< Maps tag to elements

class Domain
{
public:
    // static
    static bool PARA;     ///< Parallel code ?
    static bool WithInfo; ///< Print information ?

    // enum
    enum BryTagType_t { None_t, Element_t, Border_t, Node_t }; ///< type of boundary condition tag

    // typedefs
    typedef std::map<int,PtBCMult> MDatabase_t; ///< Map tag to M function pointer
    typedef std::map<int,Model*>   Models_t;    ///< Map tag to model pointer

    // Constructor
    Domain (Mesh::Generic const & Msh,        ///< The mesh
            Dict const          & Prps,       ///< Element properties
            Dict const          & Mdls,       ///< Model names and parameters
            Dict const          & Inis,       ///< Initial values
            char const          * FNKey=NULL, ///< Filename key to be used during the output of results
            Array<int>    const * OutV=NULL,  ///< IDs or Tags of vertices to generate output
            Array<int>    const * OutC=NULL); ///< IDs or Tags of cells to generate output

    // Destructor
    ~Domain ();

    // Methods
    void SetBCs         (Dict const & BCs);                                                                                ///< Set boundary conditions
    void PrintResults   (char const * NF="%15.6e", bool OnlySummary=false, bool WithElems=true, double Tol=1.0e-10) const; ///< Print results (Tol:tolerance to ignore zeros)
    void WriteMPY       (char const * FileKey, double SFCoef=1.0, bool PNG=false, char const * Extra=NULL) const;          ///< SFCoef: Scale-factor coefficient
    void WriteVTU       (char const * FileKey, bool DoExtrapolation=true) const;                                           ///< Write file for ParaView
    bool CheckErrorNods (Table const & NodSol, SDPair const & NodTol) const;                                               ///< Check error at nodes
    bool CheckErrorEles (Table const & EleSol, SDPair const & EleTol) const;                                               ///< Check error at the centre of elements
    bool CheckErrorIPs  (Table const & EleSol, SDPair const & EleTol) const;                                               ///< Check error at integration points of elements

    // Internal methods
    void NodalResults  () const;                                                                              ///< Calculate extrapolated values at elments IPs' to nodes
    void OutResults    (size_t IdxOut, double Time, char const * VTUFName=NULL, bool HeaderOnly=false) const; ///< Do output results
    void CalcReactions (Table & NodesReactions, SDPair & SumReactions) const;                                 ///< Calculate reactions

    // Data
    Dict          const & Prps;        ///< Element properties
    Dict          const & Inis;        ///< Initial values
    int                   NDim;        ///< Space dimension
    Models_t              Mdls;        ///< Models
    Array<Node*>          Nods;        ///< (Allocated memory) Nodes
    Array<Element*>       Eles;        ///< (Allocated memory) Elements
    Array<Node*>          TgdNods;     ///< Tagged Nodes (at boundaries)
    Array<Element*>       TgdEles;     ///< Tagged edges or faces of Elements (at boundaries)
    Array<Node*>          OutNods;     ///< Nodes for which output (results) is generated
    Array<Element*>       OutEles;     ///< Elements for which output (results) is generated
    Array<std::ofstream*> FilNods;     ///< Files with results at selected nodes (OutNods)
    Array<std::ofstream*> FilEles;     ///< Files with results at selected elements (OutEles)
    Array<Element*>       Beams;       ///< Subset of elements of type Beam
    MDatabase_t           MFuncs;      ///< Database of pointers to M functions
    Array<String>         DisplKeys;   ///< Displacement keys
    InclSupport_t         InclSupport; ///< Inclined support
    Array<Node*>          InterNodes;  ///< Nodes on the inferface between partitions (if PARA==true)
    Tag2Eles_t            Tag2Eles;    ///< Map tag to elements. Useful to activate/deactivate layers
    Array<Node*>          ActNods;     ///< Subset of active nodes (set once per call to SetBCs)
    Array<Element*>       ActEles;     ///< Subset of active elements (set once per call to SetBCs)
    Array<Node*>          NodsWithPF;  ///< Subset of nodes with prescribed F
    Array<Node*>          NodsWithPU;  ///< Subset of nodes with prescribed U
    Array<Element*>       ElesWithBC;  ///< Subset of elements with prescribed boundary conditions (ex.: beams, flowelem with convection, etc.)

    // Nodal results
    bool          HasDisps;    ///< Has displacements (ux, ...) ?
    bool          HasVeloc;    ///< Has velocities (vx, ...) ?
    Array<String> AllUKeys;    ///< All U node keys (ux, uy, H, ...)
    Array<String> AllFKeys;    ///< All F node keys (fx, fy, Q, ...)
    Array<String> AllEKeys;    ///< All element keys (sx, sy, vx, vy, ...)
    mutable Res_t NodResults;  ///< Extrapolated nodal results
    mutable Res_t NodResCount; ///< Count how many times a variable (key) was added to a node
};

bool Domain::PARA     = false;
bool Domain::WithInfo = true;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


std::ostream & operator<< (std::ostream & os, Domain const & D)
{
    String buf;
    buf.Printf("\n%s--- Models -------------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (Domain::Models_t::const_iterator p=D.Mdls.begin(); p!=D.Mdls.end(); ++p)
        os << p->first << " " << (*p->second) << std::endl;

    buf.Printf("\n%s--- Elements properties ------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    os << D.Prps << std::endl;

    buf.Printf("\n%s--- Initial values -----------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    os << D.Inis << std::endl;

    buf.Printf("\n%s--- Nodes --------------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (size_t i=0; i<D.Nods.Size(); ++i) os << (*D.Nods[i]) << "\n";

    buf.Printf("\n%s--- Elements -----------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (size_t i=0; i<D.Eles.Size(); ++i) os << (*D.Eles[i]) << "\n";

    buf.Printf("\n%s--- Elements by tags ---------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (Tag2Eles_t::const_iterator it=D.Tag2Eles.begin(); it!=D.Tag2Eles.end(); ++it)
    {
        os << "  Tag = " << it->first << " : Elements = ";
        for (size_t i=0; i<it->second.Size(); ++i) os << it->second[i]->Cell.ID << " ";
        os << "\n";
    }
    return os;
}

// Constructor and destructor

inline Domain::Domain (Mesh::Generic const & Msh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis, char const * FNKey, Array<int> const * OutV, Array<int> const * OutC)
    : Prps(ThePrps), Inis(TheInis), NDim(Msh.NDim)
{
    // info
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Domain --- allocating nodes and elements ---------------------------------------%s\n",TERM_CLR1,TERM_RST);

    // allocate models
    for (size_t i=0; i<TheMdls.Keys.Size(); ++i)
    {
        int tag = TheMdls.Keys[i];
        if (!Prps.HasKey(tag)) throw new Fatal("Domain::Domain: Prps and Mdls dictionaries must have the same tags");
        if (TheMdls(tag).HasKey("name"))
        {
            String model_name;
            MODEL.Val2Key (TheMdls(tag)("name"), model_name);
            Mdls[tag] = AllocModel (model_name, NDim, TheMdls(tag));
        }
        else throw new Fatal("Domain::Domain: Dictionary of models must have keyword 'name' defining the name of the model");
    }

    // set nodes from mesh
    std::map<int,Node*> VertID2Node;
    for (size_t i=0; i<Msh.Verts.Size(); ++i)
    {
#ifdef HAS_MPI
        if (PARA)
        {
            // skip vertices that aren't in this partition
            if (Msh.Verts[i]->PartIDs.Find(MPI::COMM_WORLD.Get_rank())<0) continue;
        }
#endif
        // new node
        Nods.Push (new Node((*Msh.Verts[i])));
        if (Msh.Verts[i]->Tag<0) TgdNods.Push (Nods.Last());
        VertID2Node[Msh.Verts[i]->ID] = Nods.Last();

        // add DOFs to node
        for (size_t k=0; k<Msh.Verts[i]->Shares.Size(); ++k)
        {
            std::pair<String,String> const & varkeys = CellTag2VarKeys (Prps, NDim, Msh.Verts[i]->Shares[k].C->Tag);
            Nods.Last()->AddDOF (varkeys.first.CStr(), varkeys.second.CStr());
        }

#ifdef HAS_MPI
        if (PARA)
        {
            // set subset of shared nodes (between partitions)
            if (Msh.Verts[i]->PartIDs.Size()>1) InterNodes.Push (Nods.Last());
        }
#endif
        // set list of nodes for output
        if (FNKey!=NULL && OutV!=NULL)
        {
            if (OutV->Find(Msh.Verts[i]->ID)>=0 || OutV->Find(Msh.Verts[i]->Tag)>=0)
            {
                String buf; buf.Printf("%s_nod_%d_%d.res", FNKey, Msh.Verts[i]->ID, Msh.Verts[i]->Tag);
                std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
                OutNods.Push (Nods.Last());
                FilNods.Push (of);
            }
        }
    }

    // check size of Nods
    if (Nods.Size()==0) throw new Fatal("FEM::Domain::Domain: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   The array of Nodes is empty.");

    // set elements from mesh
    for (size_t i=0; i<Msh.Cells.Size(); ++i)
    {
#ifdef HAS_MPI
        if (PARA)
        {
            // skip elements that aren't in this partition
            if (Msh.Cells[i]->PartID!=MPI::COMM_WORLD.Get_rank()) continue;
        }
#endif
        int tag = Msh.Cells[i]->Tag;
        if (!Prps.HasKey(tag)) throw new Fatal("Domain::SetMesh: Dictionary of element properties does not have tag=%d",tag);
        if (Prps(tag).HasKey("prob"))
        {
            // problem name
            String prob_name;
            PROB.Val2Key (Prps(tag)("prob"), prob_name);

            // model
            Model const * mdl = NULL;
            Models_t::const_iterator m = Mdls.find(tag);
            if (m!=Mdls.end()) mdl = m->second;

            // connectivity
            Array<Node*> nodes(Msh.Cells[i]->V.Size());
            for (size_t k=0; k<nodes.Size(); ++k)
            {
                std::map<int,Node*>::const_iterator it = VertID2Node.find(Msh.Cells[i]->V[k]->ID);
                if (it==VertID2Node.end()) throw new Fatal("FEM::Domain::Domain:: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   Could not find Node # %d necessary for the definition of Element # %d.",Msh.Cells[i]->V[k]->ID,Msh.Cells[i]->ID);
                nodes[k] = it->second;
            }

            // allocate element
            if (Inis.HasKey(tag)) Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), Inis(tag), nodes));
            else                  Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), SDPair() , nodes));

            // set array of Beams
            if (prob_name=="Beam") Beams.Push (Eles.Last());

            // elements with tagged borders (edges/faces)
            if (Msh.Cells[i]->BryTags.size()>0 || prob_name=="Beam") TgdEles.Push (Eles.Last());

            // map tag to element
            Tag2Eles[tag].Push (Eles.Last());

            // set list of elements for output
            if (FNKey!=NULL && OutC!=NULL)
            {
                if (OutC->Find(Msh.Cells[i]->ID)>=0 || OutC->Find(Msh.Cells[i]->Tag)>=0)
                {
                    String buf; buf.Printf("%s_ele_%d_%d.res", FNKey, Msh.Cells[i]->ID, Msh.Cells[i]->Tag);
                    std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
                    OutEles.Push (Eles.Last());
                    FilEles.Push (of);
                }
            }
        }
        else throw new Fatal("Domain::SetMesh: Dictionary of properties must have keyword 'prob' defining the type of element corresponding to a specific problem");
    }

    // check size of Eles
    if (Eles.Size()==0) throw new Fatal("FEM::Domain::Domain: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   The array of Elements is empty.");

    // displacement keys
    DisplKeys.Resize (NDim);
    DisplKeys[0] = "ux";
    DisplKeys[1] = "uy";  if (NDim==3)
    DisplKeys[2] = "uz";

    // check all available nodal data
    HasDisps = false;
    AllUKeys.Resize (0);
    AllFKeys.Resize (0);
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        for (size_t j=0; j<Nods[i]->NDOF(); ++j)
        {
            if (Nods[i]->UKey(j)=="ux") HasDisps = true;
            AllUKeys.XPush (Nods[i]->UKey(j)); // push only if item is not in array yet
            AllFKeys.XPush (Nods[i]->FKey(j)); // push only if item is not in array yet
        }
    }

    // check all available element data
    HasVeloc = false;
    Array<String> keys;
    AllEKeys.Resize (0);
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        Eles[i]->StateKeys (keys);
        for (size_t j=0; j<keys.Size(); ++j)
        {
            if (keys[j]=="vx") HasVeloc = true;
            AllEKeys.XPush (keys[j]); // push only if item is not in array yet
        }
    }

    // set maps of results
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        NodResults [Nods[i]].Resize (AllEKeys.Size()); // results at nodes
        NodResCount[Nods[i]].Resize (AllEKeys.Size()); // count how many times a varable was added to a node
    }

    // print header of output files
    OutResults (/*idxout*/0, /*time*/0, /*vtufname*/NULL, /*headeronly*/true);
}

inline Domain::~Domain()
{
    for (Domain::Models_t::iterator p=Mdls.begin(); p!=Mdls.end(); ++p) delete p->second;
    for (size_t i=0; i<Nods   .Size(); ++i) if (Nods   [i]!=NULL) delete Nods   [i];
    for (size_t i=0; i<Eles   .Size(); ++i) if (Eles   [i]!=NULL) delete Eles   [i];
    for (size_t i=0; i<FilNods.Size(); ++i) if (FilNods[i]!=NULL) { FilNods[i]->close(); delete FilNods[i]; }
    for (size_t i=0; i<FilEles.Size(); ++i) if (FilEles[i]!=NULL) { FilEles[i]->close(); delete FilEles[i]; }
}

// Methods

inline void Domain::SetBCs (Dict const & BCs)
{
    // info
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Domain --- setting BCs and activating/deactivating elements --------------------%s\n",TERM_CLR1,TERM_RST);

    // clear previous BCs
    for (size_t i=0; i<NodsWithPU.Size(); ++i) NodsWithPU[i]->DelPUs ();
    for (size_t i=0; i<NodsWithPF.Size(); ++i) NodsWithPF[i]->DelPFs ();
    for (size_t i=0; i<ElesWithBC.Size(); ++i) ElesWithBC[i]->ClrBCs ();
    NodsWithPU.Resize (0);
    NodsWithPF.Resize (0);
    ElesWithBC.Resize (0);
    InclSupport.clear ();

    // map used to track whether BC was already set or not
    std::map<int,BryTagType_t> bc_set; // maps: btag => type of BC: none, element, border, node
    for (Dict_t::const_iterator it=BCs.begin(); it!=BCs.end(); ++it) bc_set[it->first] = None_t;

    // 'BCs' at the element itself (ex: 'activate', 's':source, 'qn' at beams, etc.).
    // This must come before the setting up of 'real' boundary conditions, since, for example,
    // some nodes may become activated or deactivated
    for (Dict_t::const_iterator dict_it=BCs.begin(); dict_it!=BCs.end(); ++dict_it)
    {
        // auxiliar variables
        int            btag = dict_it->first;  // the element/boundary tag
        SDPair const & bcs  = dict_it->second; // the boundary conditions

        // check if tag corresponds to element tag (and not boundary tag)
        Tag2Eles_t::const_iterator it = Tag2Eles.find(btag);
        if (it!=Tag2Eles.end()) // found elements with this btag
        {
            // flag BC set
            bc_set[btag] = Element_t;

            // multiplier
            PtBCMult calcm = &BCMultiplier;
            if (bcs.HasKey("mfunc")) // callback specified
            {
                MDatabase_t::const_iterator q = MFuncs.find(btag);
                if (q!=MFuncs.end()) calcm = q->second;
                else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",btag);
            }

            // set BCs: ex: 'activate', 's', 'qn', etc.
            for (size_t i=0; i<it->second.Size(); ++i) // for each element
            {
                // check if bc key is not U and is not F (it cannot be, for example: 'ux', 'fx', etc., since these are applied to the edges/faces of the elements only)
                Element * ele = it->second[i];
                for (size_t k=0; k<bcs.Keys.Size(); ++k)
                {
                    if (AllUKeys.Has(bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), btag, ele->Cell.ID, ele->Cell.Tag);
                    if (AllFKeys.Has(bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), btag, ele->Cell.ID, ele->Cell.Tag);
                }

                // set BCs
                ele->SetBCs (/*idx_side(ignored)*/0, bcs, calcm);
            }
        }
    }

    // BCs at the edges/faces of elements or at nodes
    typedef std::pair<Element*,int> eleside_t;          // (element,side) pair
    typedef std::pair<int,SDPair>   tagbcs_t;           // (tag,bcs data) pair
    std::map<eleside_t,tagbcs_t>    eleside2tag_to_ubc; // maps: (ele,side) ==> U bry cond: (tag,Ubcs)
    std::map<eleside_t,tagbcs_t>    eleside2tag_to_fbc; // maps: (ele,side) ==> F bry cond: (tag,Fbcs)
    std::map<Node*,tagbcs_t>        nod2tag_to_ubc;     // map: node ==> U bry cond: (tag,Ubcs)
    std::map<Node*,tagbcs_t>        nod2tag_to_fbc;     // map: node ==> F bry cond: (tag,Fbcs)
    for (Dict_t::const_iterator dict_it=BCs.begin(); dict_it!=BCs.end(); ++dict_it)
    {
        // auxiliar variables
        int            btag  = dict_it->first;  // the element/boundary tag
        SDPair const & bcs   = dict_it->second; // the boundary conditions

        // find elements with edges equal to this btag
        for (size_t j=0; j<TgdEles.Size(); ++j) // for each element with tagged border
        {
            if (TgdEles[j]->Active) // active element
            {
                Mesh::BryTag_t const & eftags = TgdEles[j]->Cell.BryTags;
                for (Mesh::BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    int idx_side = p->first;
                    int side_tag = p->second;
                    if (side_tag==btag) // found
                    {
                        // flag BC set
                        std::map<int,BryTagType_t>::const_iterator it = bc_set.find(btag);
                        if (it->second==None_t) bc_set[btag] = Border_t;
                        else if (it->second!=Border_t) throw new Fatal("FEM::Domain::SetBCs: Boundary tag==%d cannot be applied to %s and borders (edges/faces) at the same time",btag,(it->second==Element_t?"elements":"nodes"));

                        // split Ubcs from Fbcs
                        eleside_t es(TgdEles[j],idx_side); // (edge/face,side) pair
                        for (size_t k=0; k<bcs.Keys.Size(); ++k)
                        {
                            if (AllUKeys.Has(bcs.Keys[k])) // is U key
                            {
                                eleside2tag_to_ubc[es].first = btag;
                                eleside2tag_to_ubc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else // is F key or 'qn', 'qt', etc.
                            {
                                eleside2tag_to_fbc[es].first = btag;
                                eleside2tag_to_fbc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                        }
                    }
                }
            }
        }

        // find nodes with tags equal to this btag
        for (size_t j=0; j<TgdNods.Size(); ++j)
        {
            if (TgdNods[j]->NShares>0) // active node
            {
                if (TgdNods[j]->Vert.Tag==btag) // found
                {
                    // flag BC set
                    std::map<int,BryTagType_t>::const_iterator it = bc_set.find(btag);
                    if (it->second==None_t) bc_set[btag] = Node_t;
                    else if (it->second!=Node_t) throw new Fatal("FEM::Domain::SetBCs: Boundary tag==%d cannot be applied to %s and nodes at the same time",btag,(it->second==Element_t?"elements":"borders (edges/faces)"));

                    // split U bcs from F bcs
                    for (size_t k=0; k<bcs.Keys.Size(); ++k)
                    {
                        if (bcs.Keys[k]=="inclsupport") // inclined support
                        {
                            nod2tag_to_ubc[TgdNods[j]].first  = btag;
                            nod2tag_to_ubc[TgdNods[j]].second = bcs;
                            break;
                        }
                        else
                        {
                            if (AllUKeys.Has(bcs.Keys[k])) // is U key
                            {
                                nod2tag_to_ubc[TgdNods[j]].first = btag;
                                nod2tag_to_ubc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else if (AllFKeys.Has(bcs.Keys[k])) // is F key
                            {
                                nod2tag_to_fbc[TgdNods[j]].first = btag;
                                nod2tag_to_fbc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else throw new Fatal("FEM::Domain::SetBCs: BC==%s with tag==%d cannot be specified to Node # %d", bcs.Keys[k].CStr(), btag, TgdNods[j]->Vert.ID);
                        }
                    }
                }
            }
        }
    }

    // check if all tags were set
    if (!PARA) // only in serial problems, since this processor may not know the tags of other elements/nodes in other processors
    {
        for (std::map<int,BryTagType_t>::const_iterator it=bc_set.begin(); it!=bc_set.end(); ++it)
        {
            if (it->second==None_t) throw new Fatal("FEM::Domain::SetBCs: Could not find any element, edge/face, or node with boundary tag == %d",it->first);
        }
    }

    // set F bcs at sides (edges/faces) of elements
    for (std::map<eleside_t,tagbcs_t>::iterator p=eleside2tag_to_fbc.begin(); p!=eleside2tag_to_fbc.end(); ++p)
    {
        Element      * ele      = p->first.first;
        int            idx_side = p->first.second;
        int            bc_tag   = p->second.first;
        SDPair const & bcs      = p->second.second;
        PtBCMult       calcm    = &BCMultiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        ele->SetBCs (idx_side, bcs, calcm);
    }

    // set F bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_fbc.begin(); p!=nod2tag_to_fbc.end(); ++p)
    {
        Node         * nod    = p->first;
        int            bc_tag = p->second.first;
        SDPair const & bcs    = p->second.second;
        PtBCMult       calcm  = &BCMultiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        for (StrDbl_t::const_iterator q=bcs.begin(); q!=bcs.end(); ++q)
        {
            if (q->first!="mfunc") nod->AddToPF (q->first, q->second, calcm);
        }
    }

    // set U bcs at sides (edges/faces) of elements
    for (std::map<eleside_t,tagbcs_t>::iterator p=eleside2tag_to_ubc.begin(); p!=eleside2tag_to_ubc.end(); ++p)
    {
        Element      * ele      = p->first.first;
        int            idx_side = p->first.second;
        int            bc_tag   = p->second.first;
        SDPair const & bcs      = p->second.second;
        PtBCMult       calcm    = &BCMultiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        ele->SetBCs (idx_side, bcs, calcm);
    }

    // set U bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_ubc.begin(); p!=nod2tag_to_ubc.end(); ++p)
    {
        Node         * nod    = p->first;
        int            bc_tag = p->second.first;
        SDPair const & bcs    = p->second.second;
        PtBCMult       calcm  = &BCMultiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        if (bcs.HasKey("inclsupport"))
        {
            if (NDim!=2) throw new Fatal("FEM::Domain::SetBCs: Inclined support is not implemented for 3D problems yet");
            InclSupport[nod] = bcs("alpha");
        }
        else
        {
            for (StrDbl_t::const_iterator q=bcs.begin(); q!=bcs.end(); ++q)
            {
                nod->SetPU (q->first, q->second, calcm);
            }
        }
    }

    // set subsets of active nodes and nodes with prescribed U and/or F
    ActNods.Resize (0);
    ActEles.Resize (0);
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0) ActNods   .Push (Nods[i]); 
        if (Nods[i]->NPU()>0)   NodsWithPU.Push (Nods[i]);
        if (Nods[i]->NPF()>0)   NodsWithPF.Push (Nods[i]);
    }
    for (size_t i=0; i<Eles.Size(); ++i) { if (Eles[i]->Active) ActEles.Push (Eles[i]); }
}

inline void Domain::PrintResults (char const * NF, bool OnlySummary, bool WithElems, double Tol) const
{
    printf("\n%s--- Results ------------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    // number format for text
    String nf(NF);
    size_t pos;
    pos=nf.find("g"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("g",pos+1); }
    pos=nf.find("f"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("f",pos+1); }
    pos=nf.find("e"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("e",pos+1); }

    if (!OnlySummary)
    {
        // nodes: header
        String buf;
        std::cout << TERM_CLR2 << Util::_6 << "Node";
        for (size_t i=0; i<AllUKeys.Size(); ++i) { buf.Printf(nf, AllUKeys[i].CStr());  std::cout<<buf; }
        for (size_t i=0; i<AllFKeys.Size(); ++i) { buf.Printf(nf, AllFKeys[i].CStr());  std::cout<<buf; }
        printf("%s\n",TERM_RST);

        // nodes: data
        for (size_t i=0; i<ActNods.Size(); ++i)
        {
            std::cout << Util::_6 << ActNods[i]->Vert.ID;
            for (size_t j=0; j<AllUKeys.Size(); ++j)
            {
                buf.Printf(NF, ActNods[i]->UOrZero(AllUKeys[j]));
                std::cout << buf;
            }
            for (size_t j=0; j<AllFKeys.Size(); ++j)
            {
                buf.Printf(NF, ActNods[i]->FOrZero(AllFKeys[j]));
                std::cout << buf;
            }
            printf("\n");
        }
        printf("\n");
    }

    // reactions
    String buf;
    Table  reac;
    SDPair sum;
    CalcReactions (reac, sum);
    printf("\n   Reactions\n%s %6s",TERM_CLR2,"Node");
    for (size_t i=1; i<reac.Keys.Size(); ++i)
    {
        String tmp("R");
        tmp.append (reac.Keys[i]);
        buf.Printf (nf, tmp.CStr());
        std::cout<<buf; 
    }
    printf("%s\n",TERM_RST);
    for (size_t i=0; i<reac.NRows; ++i)
    {
        printf(" %6d",static_cast<int>(reac("Node",i)));
        for (size_t j=1; j<reac.Keys.Size(); ++j) printf(NF,reac(reac.Keys[j],i));
        printf("\n");
    }
    printf("   -------------------------------------------------------------------\n");
    printf(" %6s","Sum");
    for (size_t j=0; j<sum.Keys.Size(); ++j)
    {
        double val = sum(sum.Keys[j]);
        printf(NF,(fabs(val)>Tol?val:0.0));
    }
    printf("\n");
}

inline void Domain::WriteMPY (char const * FNKey, double SFCoef, bool PNG, char const * Extra) const
{
    // bounding box
    Vec3_t min, max, del;
    min = ActNods[0]->Vert.C(0), ActNods[0]->Vert.C(1), ActNods[0]->Vert.C(2);
    max = min;
    for (size_t i=0; i<ActNods.Size(); ++i)
    {
        if (ActNods[i]->Vert.C(0)<min(0)) min(0) = ActNods[i]->Vert.C(0);
        if (ActNods[i]->Vert.C(1)<min(1)) min(1) = ActNods[i]->Vert.C(1);
        if (ActNods[i]->Vert.C(2)<min(2)) min(2) = ActNods[i]->Vert.C(2);
        if (ActNods[i]->Vert.C(0)>max(0)) max(0) = ActNods[i]->Vert.C(0);
        if (ActNods[i]->Vert.C(1)>max(1)) max(1) = ActNods[i]->Vert.C(1);
        if (ActNods[i]->Vert.C(2)>max(2)) max(2) = ActNods[i]->Vert.C(2);
    }
    del  = max - min;
    double max_dist = Norm(del);

    String fn(FNKey);  fn.append(".mpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    MPL::Header   (of);
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->Draw (of, SFCoef*max_dist);
    MPL::AddPatch (of);
    if (Extra!=NULL) of << Extra;
    if (PNG) MPL::SaveFig (FNKey, of);
    else     MPL::Show    (of);
    of.close ();
}

inline void Domain::WriteVTU (char const * FNKey, bool DoExtrapolation) const
{
    // extrapolate results
    if (DoExtrapolation) NodalResults ();

    // data
    size_t nn = ActNods.Size(); // number of active nodes
    size_t ne = ActEles.Size(); // number of active elements

    // header
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    std::map<int,int> VertID2LocID;
    for (size_t i=0; i<nn; ++i)
    {
        VertID2LocID[ActNods[i]->Vert.ID] = i;
        oss << "  " << Util::_8s <<          ActNods[i]->Vert.C(0) << " ";
        oss <<         Util::_8s <<          ActNods[i]->Vert.C(1) << " ";
        oss <<         Util::_8s << (NDim==3?ActNods[i]->Vert.C(2):0.0);
        k++;
        VTU_NEWLINE (i,k,nn,6/3-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Points>\n";

    // elements: connectivity, offsets, types
    oss << "      <Cells>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[ActEles[i]->Con.Size()] : NVertsToVTKNVerts3D[ActEles[i]->Con.Size()]);
        oss << "  ";
        for (size_t j=0; j<nne; ++j) oss << VertID2LocID[ActEles[i]->Con[j]->Vert.ID] << " ";
        k++;
        VTU_NEWLINE (i,k,ne,40/ActEles[i]->Con.Size(),oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    for (size_t i=0; i<ne; ++i)
    {
        size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[ActEles[i]->Con.Size()] : NVertsToVTKNVerts3D[ActEles[i]->Con.Size()]);
        offset += nne;
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,ne,40,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[ActEles[i]->Con.Size()];
        else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[ActEles[i]->Con.Size()];
        k++;
        VTU_NEWLINE (i,k,ne,40,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"point_scalars\">\n";
    for (size_t i=0; i<AllUKeys.Size(); ++i)
    {
        if (!(AllUKeys[i]=="ux" || AllUKeys[i]=="uy" || AllUKeys[i]=="uz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << AllUKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                oss << (k==0?"  ":" ") << ActNods[j]->UOrZero(AllUKeys[i]);
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- displacements
    if (HasDisps)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << "u" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            oss << "  " << Util::_8s <<            ActNods[i]->U("ux") << " ";
            oss <<         Util::_8s <<            ActNods[i]->U("uy") << " ";
            oss <<         Util::_8s << (NDim==3 ? ActNods[i]->U("uz") : 0.0);
            k++;
            VTU_NEWLINE (i,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- extrapolated values
    for (size_t i=0; i<AllEKeys.Size(); ++i)
    {
        if (!(AllEKeys[i]=="vx" || AllEKeys[i]=="vy" || AllEKeys[i]=="vz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << AllEKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                double cnt = NodResCount[ActNods[j]][i];
                oss << (k==0?"  ":" ") << NodResults[ActNods[j]][i] / (cnt>0.0 ? cnt : 1.0);
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- velocities
    if (HasVeloc)
    {
        size_t idx_vx=-1, idx_vy=-1, idx_vz=-1;
        for (size_t i=0; i<AllEKeys.Size(); ++i)
        {
            if (AllEKeys[i]=="vx") idx_vx = i;
            if (AllEKeys[i]=="vy") idx_vy = i;
            if (AllEKeys[i]=="vz") idx_vz = i;
        }
        oss << "        <DataArray type=\"Float32\" Name=\"" << "v" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            double cnt_x =            NodResCount[ActNods[j]][idx_vx];
            double cnt_y =            NodResCount[ActNods[j]][idx_vy];
            double cnt_z = (NDim==3 ? NodResCount[ActNods[j]][idx_vz] : 0.0);
            oss << "  " << Util::_8s <<            NodResults[ActNods[j]][idx_vx] / (cnt_x>0.0 ? cnt_x : 1.0) << " ";
            oss <<         Util::_8s <<            NodResults[ActNods[j]][idx_vy] / (cnt_y>0.0 ? cnt_y : 1.0) << " ";
            oss <<         Util::_8s << (NDim==3 ? NodResults[ActNods[j]][idx_vz] / (cnt_z>0.0 ? cnt_z : 1.0) : 0.0);
            k++;
            VTU_NEWLINE (j,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- end
    oss << "      </PointData>\n";

    // bottom
    oss << "    </Piece>\n";
    oss << "  </UnstructuredGrid>\n";
    oss << "</VTKFile>" << std::endl;

    // write to file
    String fn(FNKey); fn.append(".vtu");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline bool Domain::CheckErrorNods (Table const & NodSol, SDPair const & NodTol) const
{
    // check
    if (NodSol.NRows!=Nods.Size()) throw new Fatal("Domain::CheckErrorNods: Number of rows (%d) in table NodSol with the nodes solution must be equal to the number of nodes (%d)",NodSol.NRows,Nods.Size());

    // header
    printf("\n%s--- Error Summary --- nodes --------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

    // results
    bool error = false;

    // nodes
    Table  reac; // reactions
    SDPair sumR; // sum of reactions
    bool reactions_calculated = false;
    for (size_t i=0; i<NodSol.Keys.Size(); ++i)
    {
        // calc error
        String key = NodSol.Keys[i];
        Array<double> err(NodSol.NRows);
        if (key[0]=='R') // reaction
        {
            if (!reactions_calculated)
            {
                CalcReactions (reac, sumR);
                reactions_calculated = true;
            }
            String ukey;
            for (size_t j=1; j<key.size(); ++j) ukey.Printf("%s%c",ukey.CStr(),key[j]);
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                long pos = reac("Node").Find (Nods[j]->Vert.ID);
                if (pos<0) err[j] = 0.0;
                else       err[j] = fabs(reac(ukey,pos) - NodSol(key,j));
            }
        }
        else
        {
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                if (Nods[j]->NShares>0)
                {
                    if (AllUKeys.Has(key)) err[j] = fabs(Nods[j]->U(key.CStr()) - NodSol(key,j));
                    else                   err[j] = fabs(Nods[j]->F(key.CStr()) - NodSol(key,j));
                }
                else err[j] = 0.0;
            }
        }

        // summary
        double max_err = err.TheMax();
        double tol     = NodTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
}

inline bool Domain::CheckErrorEles (Table const & EleSol, SDPair const & EleTol) const
{
    // check
    if (EleSol.NRows!=Eles.Size()) throw new Fatal("Domain::CheckErrorEles: Number of rows (%d) in table EleSol with the elements solution must be equal to the number of elements (%d)",EleSol.NRows,Eles.Size());

    // header
    printf("\n%s--- Error Summary --- centre of elements -------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

    // results
    bool error = false;

    // elements
    for (size_t i=0; i<EleSol.Keys.Size(); ++i)
    {
        // calc error
        String key = EleSol.Keys[i];
        Array<double> err(EleSol.NRows);
        for (size_t j=0; j<Eles.Size(); ++j)
        {
            Array<SDPair> loc_res; // local results: size==number of nodes in element
            ActEles[j]->StateAtNodes (loc_res);
            SDPair ave; // average
            ave = loc_res[0];
            for (size_t k=1; k<loc_res.Size(); ++k)
            {
                ave.AddVal (key.CStr(), loc_res[k](key));
            }
            err [j] = fabs(ave(key)/loc_res.Size() - EleSol(key,j));
        }

        // summary
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
}

inline bool Domain::CheckErrorIPs (Table const & EleSol, SDPair const & EleTol) const
{
    // header
    printf("\n%s--- Error Summary --- integration points -------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

    // results
    bool error = false;

    // elements
    for (size_t i=0; i<EleSol.Keys.Size(); ++i)
    {
        // calc error
        String key = EleSol.Keys[i];
        Array<double> err(EleSol.NRows);
        for (size_t j=0; j<Eles.Size(); ++j)
        {
            if (Eles[j]->GE==NULL) throw new Fatal("Domain::CheckError: This method works only when GE (geometry element) is not NULL");
            Array<SDPair> res;
            Eles[j]->StateAtIPs (res);
            for (size_t k=0; k<Eles[j]->GE->NIP; ++k)
            {
                size_t row = k + j*Eles[j]->GE->NIP;
                err [row] = fabs(res[k](key) - EleSol(key,row));
            }
        }

        // summary
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    std::cout << "\n";

    return error;
}

// Internal methods

inline void Domain::NodalResults () const
{
    // extrapolate from elements' IPs to nodes
    for (Res_t::iterator p=NodResults .begin(); p!=NodResults .end(); ++p) p->second.SetValues(0.0);
    for (Res_t::iterator p=NodResCount.begin(); p!=NodResCount.end(); ++p) p->second.SetValues(0.0);
    for (size_t i=0; i<ActEles.Size(); ++i)
    {
        Array<SDPair> loc_res; // local results: size==number of nodes in element
        ActEles[i]->StateAtNodes (loc_res);
        for (size_t j=0; j<AllEKeys.Size(); ++j)
        {
            if (loc_res[0].HasKey(AllEKeys[j]))
            {
                for (size_t k=0; k<ActEles[i]->Con.Size(); ++k)
                {
                    NodResults [ActEles[i]->Con[k]][j] += loc_res[k](AllEKeys[j]);
                    NodResCount[ActEles[i]->Con[k]][j] += 1.0;
                }
            }
        }
    }
}

inline void Domain::OutResults (size_t IdxOut, double Time, char const * VTUFName, bool HeaderOnly) const
{
    // extrapolate values to nodes
    NodalResults ();
    
    // nodes
    for (size_t i=0; i<OutNods.Size(); ++i)
    {
        if (HeaderOnly)
        {
            (*FilNods[i]) << Util::_8s << "Time";
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllUKeys[j];
            for (size_t j=0; j<AllFKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllFKeys[j];
            for (size_t j=0; j<AllEKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllEKeys[j];
            (*FilNods[i]) << "\n";
        }
        else
        {
            (*FilNods[i]) << Util::_8s << Time;
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->UOrZero(AllUKeys[j]);
            for (size_t j=0; j<AllFKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->FOrZero(AllFKeys[j]);
            for (size_t j=0; j<AllEKeys.Size(); ++j)
            {
                double cnt = NodResCount[OutNods[i]][j];
                (*FilNods[i]) << Util::_8s << NodResults[OutNods[i]][j] / (cnt>0.0 ? cnt : 1.0);
            }
            (*FilNods[i]) << "\n";
        }
    }

    // VTU
    if (VTUFName!=NULL)
    {
        String fkey;
        fkey.Printf ("%s_%08d", VTUFName, IdxOut);
        WriteVTU    (fkey.CStr(), /*do_extrapolation*/false);
    }
}

inline void Domain::CalcReactions (Table & NodesReactions, SDPair & SumReactions) const
{
    // resize tables
    String buf("Node");
    for (size_t i=0; i<AllUKeys.Size(); ++i)
    {
        buf.Printf ("%s %s",buf.CStr(),AllUKeys[i].CStr());
        SumReactions.Set (AllUKeys[i].CStr(), 0.0);
    }
    NodesReactions.SetZero (buf.CStr(), NodsWithPU.Size());

    // find reactions
    size_t k = 0;
    for (size_t i=0; i<NodsWithPU.Size(); ++i)
    {
        Node const * nod = NodsWithPU[i];
        NodesReactions("Node", k) = nod->Vert.ID;
        std::map<String,double> R;
        nod->Reactions (R);
        for (std::map<String,double>::const_iterator it=R.begin(); it!=R.end(); ++it)
        {
            NodesReactions(it->first, k) = it->second;
            SumReactions[it->first] += it->second;
        }
        k++;
    }
}

}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
