/*****************************************************************************
 * me.c: motion estimation
 *****************************************************************************
 * Copyright (C) 2003-2020 x264 project
 *
 * Authors: Loren Merritt <lorenm@u.washington.edu>
 *          Laurent Aimar <fenrir@via.ecp.fr>
 *          Fiona Glaser <fiona@x264.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at licensing@x264.com.
 *****************************************************************************/

#include "common/common.h"
#include "macroblock.h"
#include "me.h"

/* presets selected from good points on the speed-vs-quality curve of several test videos
 * subpel_iters[i_subpel_refine] = { refine_hpel, refine_qpel, me_hpel, me_qpel }
 * where me_* are the number of EPZS iterations run on all candidate block types,
 * and refine_* are run only on the winner.
 * the subme=8,9 values are much higher because any amount of satd search makes
 * up its time by reducing the number of qpel-rd iterations. */
static const uint8_t subpel_iterations[][4] =
   {{0,0,0,0},
    {1,1,0,0},
    {0,1,1,0},
    {0,2,1,0},
    {0,2,1,1},
    {0,2,1,2},
    {0,0,2,2},
    {0,0,2,2},
    {0,0,4,10},
    {0,0,4,10},
    {0,0,4,10},
    {0,0,4,10}};

/* (x-1)%6 */
static const uint8_t mod6m1[8] = {5,0,1,2,3,4,5,0};
/* radius 2 hexagon. repeated entries are to avoid having to compute mod6 every time. */
static const int8_t hex2[8][2] = {{-1,-2}, {-2,0}, {-1,2}, {1,2}, {2,0}, {1,-2}, {-1,-2}, {-2,0}};
static const int8_t square1[9][2] = {{0,0}, {0,-1}, {0,1}, {-1,0}, {1,0}, {-1,-1}, {-1,1}, {1,-1}, {1,1}};

static void refine_subpel( x264_t *h, x264_me_t *m, int hpel_iters, int qpel_iters, int *p_halfpel_thresh, int b_refine_qpel );

#define BITS_MVD( mx, my )\
    (p_cost_mvx[(mx)<<2] + p_cost_mvy[(my)<<2])

#define COST_MV( mx, my )\
do\
{\
    int cost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE,\
                   &p_fref_w[(my)*stride+(mx)], stride )\
             + BITS_MVD(mx,my);\
    COPY3_IF_LT( bcost, cost, bmx, mx, bmy, my );\
} while( 0 )

#define COST_MV_HPEL( mx, my, cost )\
do\
{\
    intptr_t stride2 = 16;\
    pixel *src = h->mc.get_ref( pix, &stride2, m->p_fref, stride, mx, my, bw, bh, &m->weight[0] );\
    cost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE, src, stride2 )\
         + p_cost_mvx[ mx ] + p_cost_mvy[ my ];\
} while( 0 )

#define COST_MV_X3_DIR( m0x, m0y, m1x, m1y, m2x, m2y, costs )\
{\
    pixel *pix_base = p_fref_w + bmx + bmy*stride;\
    h->pixf.fpelcmp_x3[i_pixel]( p_fenc,\
        pix_base + (m0x) + (m0y)*stride,\
        pix_base + (m1x) + (m1y)*stride,\
        pix_base + (m2x) + (m2y)*stride,\
        stride, costs );\
    (costs)[0] += BITS_MVD( bmx+(m0x), bmy+(m0y) );\
    (costs)[1] += BITS_MVD( bmx+(m1x), bmy+(m1y) );\
    (costs)[2] += BITS_MVD( bmx+(m2x), bmy+(m2y) );\
}

#define COST_MV_X4_DIR( m0x, m0y, m1x, m1y, m2x, m2y, m3x, m3y, costs )\
{\
    pixel *pix_base = p_fref_w + bmx + bmy*stride;\
    h->pixf.fpelcmp_x4[i_pixel]( p_fenc,\
        pix_base + (m0x) + (m0y)*stride,\
        pix_base + (m1x) + (m1y)*stride,\
        pix_base + (m2x) + (m2y)*stride,\
        pix_base + (m3x) + (m3y)*stride,\
        stride, costs );\
    (costs)[0] += BITS_MVD( bmx+(m0x), bmy+(m0y) );\
    (costs)[1] += BITS_MVD( bmx+(m1x), bmy+(m1y) );\
    (costs)[2] += BITS_MVD( bmx+(m2x), bmy+(m2y) );\
    (costs)[3] += BITS_MVD( bmx+(m3x), bmy+(m3y) );\
}

#define COST_MV_X4( m0x, m0y, m1x, m1y, m2x, m2y, m3x, m3y )\
{\
    pixel *pix_base = p_fref_w + omx + omy*stride;\
    h->pixf.fpelcmp_x4[i_pixel]( p_fenc,\
        pix_base + (m0x) + (m0y)*stride,\
        pix_base + (m1x) + (m1y)*stride,\
        pix_base + (m2x) + (m2y)*stride,\
        pix_base + (m3x) + (m3y)*stride,\
        stride, costs );\
    costs[0] += BITS_MVD( omx+(m0x), omy+(m0y) );\
    costs[1] += BITS_MVD( omx+(m1x), omy+(m1y) );\
    costs[2] += BITS_MVD( omx+(m2x), omy+(m2y) );\
    costs[3] += BITS_MVD( omx+(m3x), omy+(m3y) );\
    COPY3_IF_LT( bcost, costs[0], bmx, omx+(m0x), bmy, omy+(m0y) );\
    COPY3_IF_LT( bcost, costs[1], bmx, omx+(m1x), bmy, omy+(m1y) );\
    COPY3_IF_LT( bcost, costs[2], bmx, omx+(m2x), bmy, omy+(m2y) );\
    COPY3_IF_LT( bcost, costs[3], bmx, omx+(m3x), bmy, omy+(m3y) );\
}

#define COST_MV_X3_ABS( m0x, m0y, m1x, m1y, m2x, m2y )\
{\
    h->pixf.fpelcmp_x3[i_pixel]( p_fenc,\
        p_fref_w + (m0x) + (m0y)*stride,\
        p_fref_w + (m1x) + (m1y)*stride,\
        p_fref_w + (m2x) + (m2y)*stride,\
        stride, costs );\
    costs[0] += p_cost_mvx[(m0x)<<2]; /* no cost_mvy */\
    costs[1] += p_cost_mvx[(m1x)<<2];\
    costs[2] += p_cost_mvx[(m2x)<<2];\
    COPY3_IF_LT( bcost, costs[0], bmx, m0x, bmy, m0y );\
    COPY3_IF_LT( bcost, costs[1], bmx, m1x, bmy, m1y );\
    COPY3_IF_LT( bcost, costs[2], bmx, m2x, bmy, m2y );\
}

/*  1  */
/* 101 */
/*  1  */
#define DIA1_ITER( mx, my )\
{\
    omx = mx; omy = my;\
    COST_MV_X4( 0,-1, 0,1, -1,0, 1,0 );\
}

#define CROSS( start, x_max, y_max )\
{\
    int i = start;\
    if( (x_max) <= X264_MIN(mv_x_max-omx, omx-mv_x_min) )\
        for( ; i < (x_max)-2; i+=4 )\
            COST_MV_X4( i,0, -i,0, i+2,0, -i-2,0 );\
    for( ; i < (x_max); i+=2 )\
    {\
        if( omx+i <= mv_x_max )\
            COST_MV( omx+i, omy );\
        if( omx-i >= mv_x_min )\
            COST_MV( omx-i, omy );\
    }\
    i = start;\
    if( (y_max) <= X264_MIN(mv_y_max-omy, omy-mv_y_min) )\
        for( ; i < (y_max)-2; i+=4 )\
            COST_MV_X4( 0,i, 0,-i, 0,i+2, 0,-i-2 );\
    for( ; i < (y_max); i+=2 )\
    {\
        if( omy+i <= mv_y_max )\
            COST_MV( omx, omy+i );\
        if( omy-i >= mv_y_min )\
            COST_MV( omx, omy-i );\
    }\
}

#define FPEL(mv) (((mv)+2)>>2) /* Convert subpel MV to fullpel with rounding... */
#define SPEL(mv) ((mv)<<2)     /* ... and the reverse. */
#define SPELx2(mv) (SPEL(mv)&0xFFFCFFFC) /* for two packed MVs */

void x264_me_search_ref( x264_t *h, x264_me_t *m, int16_t (*mvc)[2], int i_mvc, int *p_halfpel_thresh )
{
    const int bw = x264_pixel_size[m->i_pixel].w;
    const int bh = x264_pixel_size[m->i_pixel].h;
    const int i_pixel = m->i_pixel;
    const int stride = m->i_stride[0];
    int i_me_range = h->param.analyse.i_me_range;
    int bmx, bmy, bcost = COST_MAX;
    int bpred_cost = COST_MAX;
    int pmx, pmy;
    pixel *p_fenc = m->p_fenc[0];
    pixel *p_fref_w = m->p_fref_w;
    ALIGNED_ARRAY_32( pixel, pix,[16*16] );
    ALIGNED_ARRAY_8( int16_t, mvc_temp,[16],[2] );

    ALIGNED_ARRAY_16( int, costs,[16] );

    int mv_x_min = h->mb.mv_limit_fpel[0][0];
    int mv_y_min = h->mb.mv_limit_fpel[0][1];
    int mv_x_max = h->mb.mv_limit_fpel[1][0];
    int mv_y_max = h->mb.mv_limit_fpel[1][1];
/* Special version of pack to allow shortcuts in CHECK_MVRANGE */
#define pack16to32_mask2(mx,my) ((mx<<16)|(my&0x7FFF))
    uint32_t mv_min = pack16to32_mask2( -mv_x_min, -mv_y_min );
    uint32_t mv_max = pack16to32_mask2( mv_x_max, mv_y_max )|0x8000;
    uint32_t pmv, bpred_mv = 0;

#define CHECK_MVRANGE(mx,my) (!(((pack16to32_mask2(mx,my) + mv_min) | (mv_max - pack16to32_mask2(mx,my))) & 0x80004000))

    const uint16_t *p_cost_mvx = m->p_cost_mv - m->mvp[0];
    const uint16_t *p_cost_mvy = m->p_cost_mv - m->mvp[1];

    /* Try extra predictors if provided.  If subme >= 3, check subpel predictors,
     * otherwise round them to fullpel. */
    if( h->mb.i_subpel_refine >= 3 )
    {
        /* Calculate and check the MVP first */
        int bpred_mx = x264_clip3( m->mvp[0], SPEL(mv_x_min), SPEL(mv_x_max) );
        int bpred_my = x264_clip3( m->mvp[1], SPEL(mv_y_min), SPEL(mv_y_max) );
        pmv = pack16to32_mask( bpred_mx, bpred_my );
        pmx = FPEL( bpred_mx );
        pmy = FPEL( bpred_my );

        COST_MV_HPEL( bpred_mx, bpred_my, bpred_cost );
        int pmv_cost = bpred_cost;

        if( i_mvc > 0 )
        {
            /* Clip MV candidates and eliminate those equal to zero and pmv. */
            int valid_mvcs = x264_predictor_clip( mvc_temp+2, mvc, i_mvc, h->mb.mv_limit_fpel, pmv );
            if( valid_mvcs > 0 )
            {
                int i = 1, cost;
                /* We stuff pmv here to branchlessly pick between pmv and the various
                 * MV candidates. [0] gets skipped in order to maintain alignment for
                 * x264_predictor_clip. */
                M32( mvc_temp[1] ) = pmv;
                bpred_cost <<= 4;
                do
                {
                    int mx = mvc_temp[i+1][0];
                    int my = mvc_temp[i+1][1];
                    COST_MV_HPEL( mx, my, cost );
                    COPY1_IF_LT( bpred_cost, (cost << 4) + i );
                } while( ++i <= valid_mvcs );
                bpred_mx = mvc_temp[(bpred_cost&15)+1][0];
                bpred_my = mvc_temp[(bpred_cost&15)+1][1];
                bpred_cost >>= 4;
            }
        }

        /* Round the best predictor back to fullpel and get the cost, since this is where
         * we'll be starting the fullpel motion search. */
        bmx = FPEL( bpred_mx );
        bmy = FPEL( bpred_my );
        bpred_mv = pack16to32_mask(bpred_mx, bpred_my);
        if( bpred_mv&0x00030003 ) /* Only test if the tested predictor is actually subpel... */
            COST_MV( bmx, bmy );
        else                          /* Otherwise just copy the cost (we already know it) */
            bcost = bpred_cost;

        /* Test the zero vector if it hasn't been tested yet. */
        if( pmv )
        {
            if( bmx|bmy ) COST_MV( 0, 0 );
        }
        /* If a subpel mv candidate was better than the zero vector, the previous
         * fullpel check won't have gotten it even if the pmv was zero. So handle
         * that possibility here. */
        else
        {
            COPY3_IF_LT( bcost, pmv_cost, bmx, 0, bmy, 0 );
        }
    }
    else
    {
        /* Calculate and check the fullpel MVP first */
        bmx = pmx = x264_clip3( FPEL(m->mvp[0]), mv_x_min, mv_x_max );
        bmy = pmy = x264_clip3( FPEL(m->mvp[1]), mv_y_min, mv_y_max );
        pmv = pack16to32_mask( bmx, bmy );

        /* Because we are rounding the predicted motion vector to fullpel, there will be
         * an extra MV cost in 15 out of 16 cases.  However, when the predicted MV is
         * chosen as the best predictor, it is often the case that the subpel search will
         * result in a vector at or next to the predicted motion vector.  Therefore, we omit
         * the cost of the MV from the rounded MVP to avoid unfairly biasing against use of
         * the predicted motion vector.
         *
         * Disclaimer: this is a post-hoc rationalization for why this hack works. */
        bcost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE, &p_fref_w[bmy*stride+bmx], stride );

        if( i_mvc > 0 )
        {
            /* Like in subme>=3, except we also round the candidates to fullpel. */
            int valid_mvcs = x264_predictor_roundclip( mvc_temp+2, mvc, i_mvc, h->mb.mv_limit_fpel, pmv );
            if( valid_mvcs > 0 )
            {
                int i = 1, cost;
                M32( mvc_temp[1] ) = pmv;
                bcost <<= 4;
                do
                {
                    int mx = mvc_temp[i+1][0];
                    int my = mvc_temp[i+1][1];
                    cost = h->pixf.fpelcmp[i_pixel]( p_fenc, FENC_STRIDE, &p_fref_w[my*stride+mx], stride ) + BITS_MVD( mx, my );
                    COPY1_IF_LT( bcost, (cost << 4) + i );
                } while( ++i <= valid_mvcs );
                bmx = mvc_temp[(bcost&15)+1][0];
                bmy = mvc_temp[(bcost&15)+1][1];
                bcost >>= 4;
            }
        }

        /* Same as above, except the condition is simpler. */
        if( pmv )
            COST_MV( 0, 0 );
    }

    switch( h->mb.i_me_method )
    {
        case X264_ME_DIA:
        {
            /*created by yancey on 2020.5.27*/
            if(h->param.analyse._no_me_method){
                /*nothing*/
                //printf("---\n");
            }else{
                /* diamond search, radius 1 */
                bcost <<= 4;
                int i = i_me_range;
                do
                {
                    COST_MV_X4_DIR( 0,-1, 0,1, -1,0, 1,0, costs );
                    COPY1_IF_LT( bcost, (costs[0]<<4)+1 );
                    COPY1_IF_LT( bcost, (costs[1]<<4)+3 );
                    COPY1_IF_LT( bcost, (costs[2]<<4)+4 );
                    COPY1_IF_LT( bcost, (costs[3]<<4)+12 );
                    if( !(bcost&15) )
                        break;
                    bmx -= (bcost<<28)>>30;
                    bmy -= (bcost<<30)>>30;
                    bcost &= ~15;
                } while( --i && CHECK_MVRANGE(bmx, bmy) );
                bcost >>= 4;
            } 
            break;
        }
        case X264_ME_HEX:
        {
            /* hexagon search, radius 2 */
            COST_MV_X3_DIR( -2,0, -1, 2,  1, 2, costs   );
            COST_MV_X3_DIR(  2,0,  1,-2, -1,-2, costs+4 ); /* +4 for 16-byte alignment */
            bcost <<= 3;
            COPY1_IF_LT( bcost, (costs[0]<<3)+2 );
            COPY1_IF_LT( bcost, (costs[1]<<3)+3 );
            COPY1_IF_LT( bcost, (costs[2]<<3)+4 );
            COPY1_IF_LT( bcost, (costs[4]<<3)+5 );
            COPY1_IF_LT( bcost, (costs[5]<<3)+6 );
            COPY1_IF_LT( bcost, (costs[6]<<3)+7 );

            if( bcost&7 )
            {
                int dir = (bcost&7)-2;
                bmx += hex2[dir+1][0];
                bmy += hex2[dir+1][1];

                /* half hexagon, not overlapping the previous iteration */
                for( int i = (i_me_range>>1) - 1; i > 0 && CHECK_MVRANGE(bmx, bmy); i-- )
                {
                    COST_MV_X3_DIR( hex2[dir+0][0], hex2[dir+0][1],
                                    hex2[dir+1][0], hex2[dir+1][1],
                                    hex2[dir+2][0], hex2[dir+2][1],
                                    costs );
                    bcost &= ~7;
                    COPY1_IF_LT( bcost, (costs[0]<<3)+1 );
                    COPY1_IF_LT( bcost, (costs[1]<<3)+2 );
                    COPY1_IF_LT( bcost, (costs[2]<<3)+3 );
                    if( !(bcost&7) )
                        break;
                    dir += (bcost&7)-2;
                    dir = mod6m1[dir+1];
                    bmx += hex2[dir+1][0];
                    bmy += hex2[dir+1][1];
                }
            }
            bcost >>= 3;
            /* square refine */
            bcost <<= 4;
            COST_MV_X4_DIR(  0,-1,  0,1, -1,0, 1,0, costs );
            COPY1_IF_LT( bcost, (costs[0]<<4)+1 );
            COPY1_IF_LT( bcost, (costs[1]<<4)+2 );
            COPY1_IF_LT( bcost, (costs[2]<<4)+3 );
            COPY1_IF_LT( bcost, (costs[3]<<4)+4 );
            COST_MV_X4_DIR( -1,-1, -1,1, 1,-1, 1,1, costs );
            COPY1_IF_LT( bcost, (costs[0]<<4)+5 );
            COPY1_IF_LT( bcost, (costs[1]<<4)+6 );
            COPY1_IF_LT( bcost, (costs[2]<<4)+7 );
            COPY1_IF_LT( bcost, (costs[3]<<4)+8 );
            bmx += square1[bcost&15][0];
            bmy += square1[bcost&15][1];
            bcost >>= 4;
            break;
        }
    }
    /* -> qpel mv */
    uint32_t bmv = pack16to32_mask(bmx,bmy);
    uint32_t bmv_spel = SPELx2(bmv);
    if( h->mb.i_subpel_refine < 3){
        m->cost_mv = p_cost_mvx[bmx<<2] + p_cost_mvy[bmy<<2];
        m->cost = bcost;
        /* compute the real cost */
        if( bmv == pmv ) m->cost += m->cost_mv;
        M32( m->mv ) = bmv_spel;
    }
    else{
        M32(m->mv) = bpred_cost < bcost ? bpred_mv : bmv_spel;
        m->cost = X264_MIN( bpred_cost, bcost );
    }
    /* subpel refine */
    if(bcost > 666){
        if( h->mb.i_subpel_refine >= 2 ){
            int hpel = subpel_iterations[h->mb.i_subpel_refine][2];
            int qpel = subpel_iterations[h->mb.i_subpel_refine][3];
            refine_subpel( h, m, hpel, qpel, p_halfpel_thresh, 0 );
        }
    }
}
#undef COST_MV

void x264_me_refine_qpel( x264_t *h, x264_me_t *m )
{
    int hpel = subpel_iterations[h->mb.i_subpel_refine][0];
    int qpel = subpel_iterations[h->mb.i_subpel_refine][1];

    if( m->i_pixel <= PIXEL_8x8 )
        m->cost -= m->i_ref_cost;

    refine_subpel( h, m, hpel, qpel, NULL, 1 );
}

void x264_me_refine_qpel_refdupe( x264_t *h, x264_me_t *m, int *p_halfpel_thresh )
{
    refine_subpel( h, m, 0, X264_MIN( 2, subpel_iterations[h->mb.i_subpel_refine][3] ), p_halfpel_thresh, 0 );
}

#define COST_MV_SAD( mx, my ) \
{ \
    intptr_t stride = 16; \
    pixel *src = h->mc.get_ref( pix, &stride, m->p_fref, m->i_stride[0], mx, my, bw, bh, &m->weight[0] ); \
    int cost = h->pixf.fpelcmp[i_pixel]( m->p_fenc[0], FENC_STRIDE, src, stride ) \
             + p_cost_mvx[ mx ] + p_cost_mvy[ my ]; \
    COPY3_IF_LT( bcost, cost, bmx, mx, bmy, my ); \
}

#define COST_MV_SATD( mx, my, dir ) \
if( b_refine_qpel || (dir^1) != odir ) \
{ \
    intptr_t stride = 16; \
    pixel *src = h->mc.get_ref( pix, &stride, &m->p_fref[0], m->i_stride[0], mx, my, bw, bh, &m->weight[0] ); \
    int cost = h->pixf.mbcmp_unaligned[i_pixel]( m->p_fenc[0], FENC_STRIDE, src, stride ) \
             + p_cost_mvx[ mx ] + p_cost_mvy[ my ]; \
    if( b_chroma_me && cost < bcost ) \
    { \
        if( CHROMA444 ) \
        { \
            stride = 16; \
            src = h->mc.get_ref( pix, &stride, &m->p_fref[4], m->i_stride[1], mx, my, bw, bh, &m->weight[1] ); \
            cost += h->pixf.mbcmp_unaligned[i_pixel]( m->p_fenc[1], FENC_STRIDE, src, stride ); \
            if( cost < bcost ) \
            { \
                stride = 16; \
                src = h->mc.get_ref( pix, &stride, &m->p_fref[8], m->i_stride[2], mx, my, bw, bh, &m->weight[2] ); \
                cost += h->pixf.mbcmp_unaligned[i_pixel]( m->p_fenc[2], FENC_STRIDE, src, stride ); \
            } \
        } \
        else \
        { \
            h->mc.mc_chroma( pix, pix+8, 16, m->p_fref[4], m->i_stride[1], \
                             mx, 2*(my+mvy_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift ); \
            if( m->weight[1].weightfn ) \
                m->weight[1].weightfn[bw>>3]( pix, 16, pix, 16, &m->weight[1], bh>>chroma_v_shift ); \
            cost += h->pixf.mbcmp[chromapix]( m->p_fenc[1], FENC_STRIDE, pix, 16 ); \
            if( cost < bcost ) \
            { \
                if( m->weight[2].weightfn ) \
                    m->weight[2].weightfn[bw>>3]( pix+8, 16, pix+8, 16, &m->weight[2], bh>>chroma_v_shift ); \
                cost += h->pixf.mbcmp[chromapix]( m->p_fenc[2], FENC_STRIDE, pix+8, 16 ); \
            } \
        } \
    } \
    COPY4_IF_LT( bcost, cost, bmx, mx, bmy, my, bdir, dir ); \
}

static void refine_subpel( x264_t *h, x264_me_t *m, int hpel_iters, int qpel_iters, int *p_halfpel_thresh, int b_refine_qpel )
{
    const int bw = x264_pixel_size[m->i_pixel].w;
    const int bh = x264_pixel_size[m->i_pixel].h;
    const uint16_t *p_cost_mvx = m->p_cost_mv - m->mvp[0];
    const uint16_t *p_cost_mvy = m->p_cost_mv - m->mvp[1];
    const int i_pixel = m->i_pixel;
   
    ALIGNED_ARRAY_32( pixel, pix,[64*18] ); // really 17x17x2, but round up for alignment
    ALIGNED_ARRAY_16( int, costs,[4] );

    int bmx = m->mv[0];
    int bmy = m->mv[1];
    int bcost = m->cost;
    if(bcost > 365){
        if( bmy > h->mb.mv_min_spel[1] && bmy < h->mb.mv_max_spel[1] && bmx > h->mb.mv_min_spel[0] && bmx < h->mb.mv_max_spel[0] ){
            int omx = bmx, omy = bmy;
            /* We have to use mc_luma because all strides must be the same to use fpelcmp_x4 */
            h->mc.mc_luma( pix   , 64, m->p_fref, m->i_stride[0], omx, omy-1, bw, bh, &m->weight[0] );
            h->mc.mc_luma( pix+16, 64, m->p_fref, m->i_stride[0], omx, omy+1, bw, bh, &m->weight[0] );
            h->mc.mc_luma( pix+32, 64, m->p_fref, m->i_stride[0], omx-1, omy, bw, bh, &m->weight[0] );
            h->mc.mc_luma( pix+48, 64, m->p_fref, m->i_stride[0], omx+1, omy, bw, bh, &m->weight[0] );
            h->pixf.fpelcmp_x4[i_pixel]( m->p_fenc[0], pix, pix+16, pix+32, pix+48, 64, costs );
            costs[0] += p_cost_mvx[omx  ] + p_cost_mvy[omy-1];
            costs[1] += p_cost_mvx[omx  ] + p_cost_mvy[omy+1];
            costs[2] += p_cost_mvx[omx-1] + p_cost_mvy[omy  ];
            costs[3] += p_cost_mvx[omx+1] + p_cost_mvy[omy  ];
            bcost <<= 4;
            COPY1_IF_LT( bcost, (costs[0]<<4)+1 );
            COPY1_IF_LT( bcost, (costs[1]<<4)+3 );
            COPY1_IF_LT( bcost, (costs[2]<<4)+4 );
            COPY1_IF_LT( bcost, (costs[3]<<4)+12 );
            bmx -= (bcost<<28)>>30;
            bmy -= (bcost<<30)>>30;
            bcost >>= 4;
        }
    } 
    m->cost = bcost;
    m->mv[0] = bmx;
    m->mv[1] = bmy;
    m->cost_mv = p_cost_mvx[bmx] + p_cost_mvy[bmy];
}

#define BIME_CACHE( dx, dy, list )\
{\
    x264_me_t *m = m##list;\
    int i = 4 + 3*dx + dy;\
    int mvx = bm##list##x+dx;\
    int mvy = bm##list##y+dy;\
    stride[0][list][i] = bw;\
    src[0][list][i] = h->mc.get_ref( pixy_buf[list][i], &stride[0][list][i], &m->p_fref[0],\
                                     m->i_stride[0], mvx, mvy, bw, bh, x264_weight_none );\
    if( rd )\
    {\
        if( CHROMA444 )\
        {\
            stride[1][list][i] = bw;\
            src[1][list][i] = h->mc.get_ref( pixu_buf[list][i], &stride[1][list][i], &m->p_fref[4],\
                                             m->i_stride[1], mvx, mvy, bw, bh, x264_weight_none );\
            stride[2][list][i] = bw;\
            src[2][list][i] = h->mc.get_ref( pixv_buf[list][i], &stride[2][list][i], &m->p_fref[8],\
                                             m->i_stride[2], mvx, mvy, bw, bh, x264_weight_none );\
        }\
        else if( CHROMA_FORMAT )\
            h->mc.mc_chroma( pixu_buf[list][i], pixv_buf[list][i], 8, m->p_fref[4], m->i_stride[1],\
                             mvx, 2*(mvy+mv##list##y_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift );\
    }\
}

#define SATD_THRESH(cost) (cost+(cost>>4))

/* Don't unroll the BIME_CACHE loop. I couldn't find any way to force this
 * other than making its iteration count not a compile-time constant. */
#define x264_iter_kludge x264_template(iter_kludge)
int x264_iter_kludge = 0;

static ALWAYS_INLINE void me_refine_bidir( x264_t *h, x264_me_t *m0, x264_me_t *m1, int i_weight, int i8, int i_lambda2, int rd )
{
    int x = i8&1;
    int y = i8>>1;
    int s8 = X264_SCAN8_0 + 2*x + 16*y;
    int16_t *cache0_mv = h->mb.cache.mv[0][s8];
    int16_t *cache1_mv = h->mb.cache.mv[1][s8];
    const int i_pixel = m0->i_pixel;
    const int bw = x264_pixel_size[i_pixel].w;
    const int bh = x264_pixel_size[i_pixel].h;
    ALIGNED_ARRAY_32( pixel, pixy_buf,[2],[9][16*16] );
    ALIGNED_ARRAY_32( pixel, pixu_buf,[2],[9][16*16] );
    ALIGNED_ARRAY_32( pixel, pixv_buf,[2],[9][16*16] );
    pixel *src[3][2][9];
    int chromapix = h->luma2chroma_pixel[i_pixel];
    int chroma_v_shift = CHROMA_V_SHIFT;
    int chroma_x = (8 >> CHROMA_H_SHIFT) * x;
    int chroma_y = (8 >> chroma_v_shift) * y;
    pixel *pix  = &h->mb.pic.p_fdec[0][8*x + 8*y*FDEC_STRIDE];
    pixel *pixu = CHROMA_FORMAT ? &h->mb.pic.p_fdec[1][chroma_x + chroma_y*FDEC_STRIDE] : NULL;
    pixel *pixv = CHROMA_FORMAT ? &h->mb.pic.p_fdec[2][chroma_x + chroma_y*FDEC_STRIDE] : NULL;
    int ref0 = h->mb.cache.ref[0][s8];
    int ref1 = h->mb.cache.ref[1][s8];
    const int mv0y_offset = chroma_v_shift & MB_INTERLACED & ref0 ? (h->mb.i_mb_y & 1)*4 - 2 : 0;
    const int mv1y_offset = chroma_v_shift & MB_INTERLACED & ref1 ? (h->mb.i_mb_y & 1)*4 - 2 : 0;
    intptr_t stride[3][2][9];
    int bm0x = m0->mv[0];
    int bm0y = m0->mv[1];
    int bm1x = m1->mv[0];
    int bm1y = m1->mv[1];
    int bcost = COST_MAX;
    int mc_list0 = 1, mc_list1 = 1;
    uint64_t bcostrd = COST_MAX64;
    uint16_t amvd;
    /* each byte of visited represents 8 possible m1y positions, so a 4D array isn't needed */
    ALIGNED_ARRAY_64( uint8_t, visited,[8],[8][8] );
    /* all permutations of an offset in up to 2 of the dimensions */
    ALIGNED_4( static const int8_t dia4d[33][4] ) =
    {
        {0,0,0,0},
        {0,0,0,1}, {0,0,0,-1}, {0,0,1,0}, {0,0,-1,0},
        {0,1,0,0}, {0,-1,0,0}, {1,0,0,0}, {-1,0,0,0},
        {0,0,1,1}, {0,0,-1,-1},{0,1,1,0}, {0,-1,-1,0},
        {1,1,0,0}, {-1,-1,0,0},{1,0,0,1}, {-1,0,0,-1},
        {0,1,0,1}, {0,-1,0,-1},{1,0,1,0}, {-1,0,-1,0},
        {0,0,-1,1},{0,0,1,-1}, {0,-1,1,0},{0,1,-1,0},
        {-1,1,0,0},{1,-1,0,0}, {1,0,0,-1},{-1,0,0,1},
        {0,-1,0,1},{0,1,0,-1}, {-1,0,1,0},{1,0,-1,0},
    };

    if( bm0y < h->mb.mv_min_spel[1] + 8 || bm1y < h->mb.mv_min_spel[1] + 8 ||
        bm0y > h->mb.mv_max_spel[1] - 8 || bm1y > h->mb.mv_max_spel[1] - 8 ||
        bm0x < h->mb.mv_min_spel[0] + 8 || bm1x < h->mb.mv_min_spel[0] + 8 ||
        bm0x > h->mb.mv_max_spel[0] - 8 || bm1x > h->mb.mv_max_spel[0] - 8 )
        return;

    if( rd && m0->i_pixel != PIXEL_16x16 && i8 != 0 )
    {
        x264_mb_predict_mv( h, 0, i8<<2, bw>>2, m0->mvp );
        x264_mb_predict_mv( h, 1, i8<<2, bw>>2, m1->mvp );
    }

    const uint16_t *p_cost_m0x = m0->p_cost_mv - m0->mvp[0];
    const uint16_t *p_cost_m0y = m0->p_cost_mv - m0->mvp[1];
    const uint16_t *p_cost_m1x = m1->p_cost_mv - m1->mvp[0];
    const uint16_t *p_cost_m1y = m1->p_cost_mv - m1->mvp[1];

    h->mc.memzero_aligned( visited, sizeof(uint8_t[8][8][8]) );

    for( int pass = 0; pass < 8; pass++ )
    {
        int bestj = 0;
        /* check all mv pairs that differ in at most 2 components from the current mvs. */
        /* doesn't do chroma ME. this probably doesn't matter, as the gains
         * from bidir ME are the same with and without chroma ME. */

        if( mc_list0 )
            for( int j = x264_iter_kludge; j < 9; j++ )
                BIME_CACHE( square1[j][0], square1[j][1], 0 );

        if( mc_list1 )
            for( int j = x264_iter_kludge; j < 9; j++ )
                BIME_CACHE( square1[j][0], square1[j][1], 1 );

        for( int j = !!pass; j < 33; j++ )
        {
            int m0x = dia4d[j][0] + bm0x;
            int m0y = dia4d[j][1] + bm0y;
            int m1x = dia4d[j][2] + bm1x;
            int m1y = dia4d[j][3] + bm1y;
            if( !pass || !((visited[(m0x)&7][(m0y)&7][(m1x)&7] & (1<<((m1y)&7)))) )
            {
                int i0 = 4 + 3*dia4d[j][0] + dia4d[j][1];
                int i1 = 4 + 3*dia4d[j][2] + dia4d[j][3];
                visited[(m0x)&7][(m0y)&7][(m1x)&7] |= (1<<((m1y)&7));
                h->mc.avg[i_pixel]( pix, FDEC_STRIDE, src[0][0][i0], stride[0][0][i0], src[0][1][i1], stride[0][1][i1], i_weight );
                int cost = h->pixf.mbcmp[i_pixel]( m0->p_fenc[0], FENC_STRIDE, pix, FDEC_STRIDE )
                         + p_cost_m0x[m0x] + p_cost_m0y[m0y] + p_cost_m1x[m1x] + p_cost_m1y[m1y];
                if( rd )
                {
                    if( cost < SATD_THRESH(bcost) )
                    {
                        bcost = X264_MIN( cost, bcost );
                        M32( cache0_mv ) = pack16to32_mask(m0x,m0y);
                        M32( cache1_mv ) = pack16to32_mask(m1x,m1y);
                        if( CHROMA444 )
                        {
                            h->mc.avg[i_pixel]( pixu, FDEC_STRIDE, src[1][0][i0], stride[1][0][i0], src[1][1][i1], stride[1][1][i1], i_weight );
                            h->mc.avg[i_pixel]( pixv, FDEC_STRIDE, src[2][0][i0], stride[2][0][i0], src[2][1][i1], stride[2][1][i1], i_weight );
                        }
                        else if( CHROMA_FORMAT )
                        {
                            h->mc.avg[chromapix]( pixu, FDEC_STRIDE, pixu_buf[0][i0], 8, pixu_buf[1][i1], 8, i_weight );
                            h->mc.avg[chromapix]( pixv, FDEC_STRIDE, pixv_buf[0][i0], 8, pixv_buf[1][i1], 8, i_weight );
                        }
                        uint64_t costrd = x264_rd_cost_part( h, i_lambda2, i8*4, m0->i_pixel );
                        COPY2_IF_LT( bcostrd, costrd, bestj, j );
                    }
                }
                else
                    COPY2_IF_LT( bcost, cost, bestj, j );
            }
        }

        if( !bestj )
            break;

        bm0x += dia4d[bestj][0];
        bm0y += dia4d[bestj][1];
        bm1x += dia4d[bestj][2];
        bm1y += dia4d[bestj][3];

        mc_list0 = M16( &dia4d[bestj][0] );
        mc_list1 = M16( &dia4d[bestj][2] );
    }

    if( rd )
    {
        x264_macroblock_cache_mv ( h, 2*x, 2*y, bw>>2, bh>>2, 0, pack16to32_mask(bm0x, bm0y) );
        amvd = pack8to16( X264_MIN(abs(bm0x - m0->mvp[0]),33), X264_MIN(abs(bm0y - m0->mvp[1]),33) );
        x264_macroblock_cache_mvd( h, 2*x, 2*y, bw>>2, bh>>2, 0, amvd );

        x264_macroblock_cache_mv ( h, 2*x, 2*y, bw>>2, bh>>2, 1, pack16to32_mask(bm1x, bm1y) );
        amvd = pack8to16( X264_MIN(abs(bm1x - m1->mvp[0]),33), X264_MIN(abs(bm1y - m1->mvp[1]),33) );
        x264_macroblock_cache_mvd( h, 2*x, 2*y, bw>>2, bh>>2, 1, amvd );
    }

    m0->mv[0] = bm0x;
    m0->mv[1] = bm0y;
    m1->mv[0] = bm1x;
    m1->mv[1] = bm1y;
}

void x264_me_refine_bidir_satd( x264_t *h, x264_me_t *m0, x264_me_t *m1, int i_weight )
{
    me_refine_bidir( h, m0, m1, i_weight, 0, 0, 0 );
}

void x264_me_refine_bidir_rd( x264_t *h, x264_me_t *m0, x264_me_t *m1, int i_weight, int i8, int i_lambda2 )
{
    /* Motion compensation is done as part of bidir_rd; don't repeat
     * it in encoding. */
    h->mb.b_skip_mc = 1;
    me_refine_bidir( h, m0, m1, i_weight, i8, i_lambda2, 1 );
    h->mb.b_skip_mc = 0;
}

#undef COST_MV_SATD
#define COST_MV_SATD( mx, my, dst, avoid_mvp ) \
{ \
    if( !avoid_mvp || !(mx == pmx && my == pmy) ) \
    { \
        h->mc.mc_luma( pix, FDEC_STRIDE, m->p_fref, m->i_stride[0], mx, my, bw, bh, &m->weight[0] ); \
        dst = h->pixf.mbcmp[i_pixel]( m->p_fenc[0], FENC_STRIDE, pix, FDEC_STRIDE ) \
            + p_cost_mvx[mx] + p_cost_mvy[my]; \
        COPY1_IF_LT( bsatd, dst ); \
    } \
    else \
        dst = COST_MAX; \
}

#define COST_MV_RD( mx, my, satd, do_dir, mdir ) \
{ \
    if( satd <= SATD_THRESH(bsatd) ) \
    { \
        uint64_t cost; \
        M32( cache_mv ) = pack16to32_mask(mx,my); \
        if( CHROMA444 ) \
        { \
            h->mc.mc_luma( pixu, FDEC_STRIDE, &m->p_fref[4], m->i_stride[1], mx, my, bw, bh, &m->weight[1] ); \
            h->mc.mc_luma( pixv, FDEC_STRIDE, &m->p_fref[8], m->i_stride[2], mx, my, bw, bh, &m->weight[2] ); \
        } \
        else if( CHROMA_FORMAT && m->i_pixel <= PIXEL_8x8 ) \
        { \
            h->mc.mc_chroma( pixu, pixv, FDEC_STRIDE, m->p_fref[4], m->i_stride[1], \
                             mx, 2*(my+mvy_offset)>>chroma_v_shift, bw>>1, bh>>chroma_v_shift ); \
            if( m->weight[1].weightfn ) \
                m->weight[1].weightfn[bw>>3]( pixu, FDEC_STRIDE, pixu, FDEC_STRIDE, &m->weight[1], bh>>chroma_v_shift ); \
            if( m->weight[2].weightfn ) \
                m->weight[2].weightfn[bw>>3]( pixv, FDEC_STRIDE, pixv, FDEC_STRIDE, &m->weight[2], bh>>chroma_v_shift ); \
        } \
        cost = x264_rd_cost_part( h, i_lambda2, i4, m->i_pixel ); \
        COPY4_IF_LT( bcost, cost, bmx, mx, bmy, my, dir, do_dir?mdir:dir ); \
    } \
}

void x264_me_refine_qpel_rd( x264_t *h, x264_me_t *m, int i_lambda2, int i4, int i_list )
{
    int16_t *cache_mv = h->mb.cache.mv[i_list][x264_scan8[i4]];
    const uint16_t *p_cost_mvx, *p_cost_mvy;
    const int bw = x264_pixel_size[m->i_pixel].w;
    const int bh = x264_pixel_size[m->i_pixel].h;
    const int i_pixel = m->i_pixel;
    int chroma_v_shift = CHROMA_V_SHIFT;
    int mvy_offset = chroma_v_shift & MB_INTERLACED & m->i_ref ? (h->mb.i_mb_y & 1)*4 - 2 : 0;

    uint64_t bcost = COST_MAX64;
    int bmx = m->mv[0];
    int bmy = m->mv[1];
    int omx, omy, pmx, pmy;
    int satd, bsatd;
    int dir = -2;
    int i8 = i4>>2;
    uint16_t amvd;

    pixel *pix  = &h->mb.pic.p_fdec[0][block_idx_xy_fdec[i4]];
    pixel *pixu, *pixv;
    if( CHROMA444 )
    {
        pixu = &h->mb.pic.p_fdec[1][block_idx_xy_fdec[i4]];
        pixv = &h->mb.pic.p_fdec[2][block_idx_xy_fdec[i4]];
    }
    else if( CHROMA_FORMAT )
    {
        pixu = &h->mb.pic.p_fdec[1][(i8>>1)*(8*FDEC_STRIDE>>chroma_v_shift)+(i8&1)*4];
        pixv = &h->mb.pic.p_fdec[2][(i8>>1)*(8*FDEC_STRIDE>>chroma_v_shift)+(i8&1)*4];
    }
    else
    {
        pixu = NULL;
        pixv = NULL;
    }

    h->mb.b_skip_mc = 1;

    if( m->i_pixel != PIXEL_16x16 && i4 != 0 )
        x264_mb_predict_mv( h, i_list, i4, bw>>2, m->mvp );
    pmx = m->mvp[0];
    pmy = m->mvp[1];
    p_cost_mvx = m->p_cost_mv - pmx;
    p_cost_mvy = m->p_cost_mv - pmy;
    COST_MV_SATD( bmx, bmy, bsatd, 0 );
    if( m->i_pixel != PIXEL_16x16 )
        COST_MV_RD( bmx, bmy, 0, 0, 0 )
    else
        bcost = m->cost;

    /* check the predicted mv */
    if( (bmx != pmx || bmy != pmy)
        && pmx >= h->mb.mv_min_spel[0] && pmx <= h->mb.mv_max_spel[0]
        && pmy >= h->mb.mv_min_spel[1] && pmy <= h->mb.mv_max_spel[1] )
    {
        COST_MV_SATD( pmx, pmy, satd, 0 );
        COST_MV_RD  ( pmx, pmy, satd, 0, 0 );
        /* The hex motion search is guaranteed to not repeat the center candidate,
         * so if pmv is chosen, set the "MV to avoid checking" to bmv instead. */
        if( bmx == pmx && bmy == pmy )
        {
            pmx = m->mv[0];
            pmy = m->mv[1];
        }
    }

    if( bmy < h->mb.mv_min_spel[1] + 3 || bmy > h->mb.mv_max_spel[1] - 3 ||
        bmx < h->mb.mv_min_spel[0] + 3 || bmx > h->mb.mv_max_spel[0] - 3 )
    {
        h->mb.b_skip_mc = 0;
        return;
    }

    /* subpel hex search, same pattern as ME HEX. */
    dir = -2;
    omx = bmx;
    omy = bmy;
    for( int j = 0; j < 6; j++ )
    {
        COST_MV_SATD( omx + hex2[j+1][0], omy + hex2[j+1][1], satd, 1 );
        COST_MV_RD  ( omx + hex2[j+1][0], omy + hex2[j+1][1], satd, 1, j );
    }

    if( dir != -2 )
    {
        /* half hexagon, not overlapping the previous iteration */
        for( int i = 1; i < 10; i++ )
        {
            const int odir = mod6m1[dir+1];
            if( bmy < h->mb.mv_min_spel[1] + 3 ||
                bmy > h->mb.mv_max_spel[1] - 3 )
                break;
            dir = -2;
            omx = bmx;
            omy = bmy;
            for( int j = 0; j < 3; j++ )
            {
                COST_MV_SATD( omx + hex2[odir+j][0], omy + hex2[odir+j][1], satd, 1 );
                COST_MV_RD  ( omx + hex2[odir+j][0], omy + hex2[odir+j][1], satd, 1, odir-1+j );
            }
            if( dir == -2 )
                break;
        }
    }

    /* square refine, same pattern as ME HEX. */
    omx = bmx;
    omy = bmy;
    for( int i = 0; i < 8; i++ )
    {
        COST_MV_SATD( omx + square1[i+1][0], omy + square1[i+1][1], satd, 1 );
        COST_MV_RD  ( omx + square1[i+1][0], omy + square1[i+1][1], satd, 0, 0 );
    }

    m->cost = bcost;
    m->mv[0] = bmx;
    m->mv[1] = bmy;
    x264_macroblock_cache_mv ( h, block_idx_x[i4], block_idx_y[i4], bw>>2, bh>>2, i_list, pack16to32_mask(bmx, bmy) );
    amvd = pack8to16( X264_MIN(abs(bmx - m->mvp[0]),66), X264_MIN(abs(bmy - m->mvp[1]),66) );
    x264_macroblock_cache_mvd( h, block_idx_x[i4], block_idx_y[i4], bw>>2, bh>>2, i_list, amvd );
    h->mb.b_skip_mc = 0;
}
